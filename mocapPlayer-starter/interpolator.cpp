#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "motion.h"
#include "interpolator.h"
#include "types.h"
#include "transform.h"
#include "quaternion.h"

Interpolator::Interpolator()
{
  //Set default interpolation type
  m_InterpolationType = LINEAR;

  //set default angle representation to use for interpolation
  m_AngleRepresentation = EULER;
}

Interpolator::~Interpolator()
{
}

//Create interpolated motion
void Interpolator::Interpolate(Motion * pInputMotion, Motion ** pOutputMotion, int N)
{
  //Allocate new motion
  *pOutputMotion = new Motion(pInputMotion->GetNumFrames(), pInputMotion->GetSkeleton());

  //Perform the interpolation
  if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == EULER))
    LinearInterpolationEuler(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == QUATERNION))
    LinearInterpolationQuaternion(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == EULER))
    BezierInterpolationEuler(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == QUATERNION))
    BezierInterpolationQuaternion(pInputMotion, *pOutputMotion, N);
  else
  {
    printf("Error: unknown interpolation / angle representation type.\n");
    exit(1);
  }
}

void Interpolator::Rotation2Euler(double R[9], double angles[3])
{
  double cy = sqrt(R[0]*R[0] + R[3]*R[3]);

  if (cy > 16*DBL_EPSILON)
  {
    angles[0] = atan2(R[7], R[8]);
    angles[1] = atan2(-R[6], cy);
    angles[2] = atan2(R[3], R[0]);
  }
  else
  {
    angles[0] = atan2(-R[5], R[4]);
    angles[1] = atan2(-R[6], cy);
    angles[2] = 0;
  }

  for(int i=0; i<3; i++)
    angles[i] *= 180 / M_PI;
}

void Interpolator::Euler2Rotation(double angles[3], double R[9])
{
  //converting angles into radians
  double angles_radians[3];
  for(int i=0; i<3; i++)
    angles_radians[i] = angles[i] * M_PI / 180;

  double R1[3][3], R2[3][3], R3[3][3], temp[3][3], finalR[3][3];

  R3[0][0] = cos(angles_radians[2]);       R3[0][1] = -1 * sin(angles_radians[2]);  R3[0][2] = 0;
  R3[1][0] = sin(angles_radians[2]);       R3[1][1] = cos(angles_radians[2]);       R3[1][2] = 0;
  R3[2][0] = 0;                            R3[2][1] = 0;                            R3[2][2] = 1;

  R2[0][0] = cos(angles_radians[1]);       R2[0][1] = 0;                            R2[0][2] = sin(angles_radians[1]);
  R2[1][0] = 0;                            R2[1][1] = 1;                            R2[1][2] = 0;
  R2[2][0] =  -1 * sin(angles_radians[1]); R2[2][1] = 0;                            R2[2][2] = cos(angles_radians[1]);

  R1[0][0] = 1;                            R1[0][1] = 0;                            R1[0][2] = 0;
  R1[1][0] = 0;                            R1[1][1] = cos(angles_radians[0]);       R1[1][2] = -1 * sin(angles_radians[0]);
  R1[2][0] = 0;                            R1[2][1] = sin(angles_radians[0]);       R1[2][2] = cos(angles_radians[0]);

  //R = R3 * R2 * R1;
  matrix_mult_3x3(R3, R2, temp);
  matrix_mult_3x3(temp, R1, finalR);

  R[0] = finalR[0][0];
  R[1] = finalR[0][1];
  R[2] = finalR[0][2];
  R[3] = finalR[1][0];
  R[4] = finalR[1][1];
  R[5] = finalR[1][2];
  R[6] = finalR[2][0];
  R[7] = finalR[2][1];
  R[8] = finalR[2][2];
}

void Interpolator::Euler2Quaternion(double angles[3], Quaternion<double> & q)
{
  //converting Euler angles to Rotation Matrix and then to Quaternion
  double rotMatrix[9];
  Euler2Rotation(angles, rotMatrix);
  q = Quaternion<double>::Matrix2Quaternion(rotMatrix);
  q.Normalize();
}

void Interpolator::Quaternion2Euler(Quaternion<double> & q, double angles[3])
{
  //converting Quaternion to Rotation Matrix and then to Euler angles
  double rotMatrix[9];
  q.Normalize();
  q.Quaternion2Matrix(rotMatrix);
  Rotation2Euler(rotMatrix, angles);
}

vector Interpolator::EulerSlerp(double t, vector pStart, vector pEnd)
{
  vector result;
  result = pStart * (1 - t) + pEnd * t;
  return result;
}

vector Interpolator::EulerDouble(vector p, vector q)
{
  vector result;
  result = EulerSlerp(2.0, p, q);
  return result;
}

Quaternion<double> Interpolator::Slerp(double t, Quaternion<double> qStart, Quaternion<double> qEnd)
{
  Quaternion<double> result;

  double cos_theta = qStart.Gets() * qEnd.Gets() + qStart.Getx() * qEnd.Getx() +
                     qStart.Gety() * qEnd.Gety() + qStart.Getz() * qEnd.Getz();

  double theta, c1, c2;

  if (cos_theta < 0)
  {
    cos_theta *= -1;
    qEnd = qEnd * (-1.0);
  }
  if (fabs(cos_theta - 1) < 1E-10)
    theta = 0;
  else
    theta = acos(cos_theta);

  //theta = 0 means that there is no rotation, so we return qStart
  if (sin(theta) == 0)
  {
    return qStart;
  }

  if (fabs(sin(theta)) < 1E-10)
  {
    result = (1 - t) * qStart + t * qEnd;
  }
  else
  {
    c1 = sin((1 - t) * theta) / sin(theta);
    c2 = sin(t * theta) / sin(theta);
    result = (c1 * qStart) + (c2 * qEnd);
  }
  result.Normalize();
  return result;
}

Quaternion<double> Interpolator::Double(Quaternion<double> p, Quaternion<double> q)
{
  Quaternion<double> result;
  result = Slerp(2.0, p, q);
  return result;
}

vector Interpolator::DeCasteljauEuler(double t, vector p0, vector p1, vector p2, vector p3)
{
  vector result;
  vector q0, q1, q2, r0, r1;

  q0 = EulerSlerp(t, p0, p1);
  q1 = EulerSlerp(t, p1, p2);
  q2 = EulerSlerp(t, p2, p3);

  r0 = EulerSlerp(t, q0, q1);
  r1 = EulerSlerp(t, q1, q2);

  result = EulerSlerp(t, r0, r1);

  return result;
}

Quaternion<double> Interpolator::DeCasteljauQuaternion(double t, Quaternion<double> p0, Quaternion<double> p1, Quaternion<double> p2, Quaternion<double> p3)
{
  Quaternion<double> result;
  Quaternion<double> q0, q1, q2, r0, r1;

  q0 = Slerp(t, p0, p1);
  q1 = Slerp(t, p1, p2);
  q2 = Slerp(t, p2, p3);

  r0 = Slerp(t, q0, q1);
  r1 = Slerp(t, q1, q2);

  result = Slerp(t, r0, r1);

  return result;
}

void Interpolator::LinearInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

  int startKeyframe = 0;
  while (startKeyframe + N + 1 < inputLength)
  {
    int endKeyframe = startKeyframe + N + 1;

    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    // interpolate in between
    for(int frame=1; frame<=N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N+1);

      // interpolate root position
      interpolatedPosture.root_pos = startPosture->root_pos * (1-t) + endPosture->root_pos * t;

      // interpolate bone rotations
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
        interpolatedPosture.bone_rotation[bone] = startPosture->bone_rotation[bone] * (1-t) + endPosture->bone_rotation[bone] * t;

      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;
  }

  for(int frame=startKeyframe+1; frame<inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::BezierInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1
  vector a, b, an_;
  vector pStart, pEnd, pA, pB;

  int startKeyframe = 0;
  while (startKeyframe + N + 1 < inputLength)
  {
    int endKeyframe = startKeyframe + N + 1;
    Posture * prevPosture = NULL;
    Posture * nextPosture = NULL;

    if (startKeyframe - N - 1 >= 0)
    {
      prevPosture = pInputMotion->GetPosture(startKeyframe - N - 1);
    }
    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);
    if (endKeyframe + N + 1 < inputLength)
    {
      nextPosture = pInputMotion->GetPosture(endKeyframe + N + 1);
    }

    if (prevPosture == NULL)
    {
      pA = EulerSlerp(1.0/3, startPosture->root_pos, EulerDouble(nextPosture->root_pos, endPosture->root_pos));

      an_ = EulerSlerp(0.5, EulerDouble(startPosture->root_pos, endPosture->root_pos), nextPosture->root_pos);
      pB = EulerSlerp(-1.0/3, endPosture->root_pos, an_);
    }
    else if (nextPosture == NULL)
    {
      an_ = EulerSlerp(1.0/2, EulerDouble(prevPosture->root_pos, startPosture->root_pos), endPosture->root_pos);
      pA = EulerSlerp(1.0/3, startPosture->root_pos, an_);

      pB = EulerSlerp(1.0/3, endPosture->root_pos, EulerDouble(prevPosture->root_pos, startPosture->root_pos));
    }
    else
    {
      an_ = EulerSlerp(1.0/2, EulerDouble(prevPosture->root_pos, startPosture->root_pos), endPosture->root_pos);
      pA = EulerSlerp(1.0/3, startPosture->root_pos, an_);

      an_ = EulerSlerp(1.0/2, EulerDouble(startPosture->root_pos, endPosture->root_pos), nextPosture->root_pos);
      pB = EulerSlerp(-1.0/3, endPosture->root_pos, an_);
    }

    pStart = startPosture->root_pos;
    pEnd = endPosture->root_pos;

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    // interpolate in between
    for(int frame=1; frame<=N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N+1);

      // interpolate root position
      interpolatedPosture.root_pos = DeCasteljauEuler(t, pStart, pA, pB, pEnd);

      // interpolate bone rotations
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
      {
        if (prevPosture == NULL)
        {
          a = EulerSlerp(1.0/3, startPosture->bone_rotation[bone], EulerDouble(nextPosture->bone_rotation[bone], endPosture->bone_rotation[bone]));

          an_ = EulerSlerp(0.5, EulerDouble(startPosture->bone_rotation[bone], endPosture->bone_rotation[bone]), nextPosture->bone_rotation[bone]);
          b = EulerSlerp(-1.0/3, endPosture->bone_rotation[bone], an_);
        }
        else if (nextPosture == NULL)
        {
          an_ = EulerSlerp(1.0/2, EulerDouble(prevPosture->bone_rotation[bone], startPosture->bone_rotation[bone]), endPosture->bone_rotation[bone]);
          a = EulerSlerp(1.0/3, endPosture->bone_rotation[bone], an_);

          b = EulerSlerp(1.0/3, endPosture->bone_rotation[bone], EulerDouble(prevPosture->bone_rotation[bone], startPosture->bone_rotation[bone]));
        }
        else
        {
          an_ = EulerSlerp(1.0/2, EulerDouble(prevPosture->bone_rotation[bone], startPosture->bone_rotation[bone]), endPosture->bone_rotation[bone]);
          a = EulerSlerp(1.0/3, startPosture->bone_rotation[bone], an_);

          an_ = EulerSlerp(1.0/2, EulerDouble(startPosture->bone_rotation[bone], endPosture->bone_rotation[bone]), nextPosture->bone_rotation[bone]);
          b = EulerSlerp(-1.0/3, endPosture->bone_rotation[bone], an_);
        }

        interpolatedPosture.bone_rotation[bone] = DeCasteljauEuler(t, startPosture->bone_rotation[bone], a, b, endPosture->bone_rotation[bone]);
      }
      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;
  }

  for(int frame=startKeyframe+1; frame<inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::LinearInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

  int startKeyframe = 0;

  while (startKeyframe + N + 1 < inputLength)
  {
    int endKeyframe = startKeyframe + N + 1;

    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    // interpolate in between
    for(int frame=1; frame<=N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N+1);
      Quaternion<double> interpolatedPostureQuaternion, pStart, pEnd;

      // interpolate root position linearly
      interpolatedPosture.root_pos = startPosture->root_pos * (1-t) + endPosture->root_pos * t;

      // interpolate bone rotations
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
      {
        Euler2Quaternion(startPosture->bone_rotation[bone].p, pStart);
        Euler2Quaternion(endPosture->bone_rotation[bone].p, pEnd);
        interpolatedPostureQuaternion = Slerp(t, pStart, pEnd);
        Quaternion2Euler(interpolatedPostureQuaternion, interpolatedPosture.bone_rotation[bone].p);
      }

      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;
  }

  for(int frame=startKeyframe+1; frame<inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::BezierInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1
  vector an_;
  vector pStart, pEnd, pA, pB;

  int startKeyframe = 0;
  while (startKeyframe + N + 1 < inputLength)
  {
    int endKeyframe = startKeyframe + N + 1;
    Posture * prevPosture = NULL;
    Posture * nextPosture = NULL;

    if (startKeyframe - N - 1 >= 0)
    {
      prevPosture = pInputMotion->GetPosture(startKeyframe - N - 1);
    }
    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);
    if (endKeyframe + N + 1 < inputLength)
    {
      nextPosture = pInputMotion->GetPosture(endKeyframe + N + 1);
    }

    if (prevPosture == NULL)
    {
      pA = EulerSlerp(1.0/3, startPosture->root_pos, EulerDouble(nextPosture->root_pos, endPosture->root_pos));

      an_ = EulerSlerp(0.5, EulerDouble(startPosture->root_pos, endPosture->root_pos), nextPosture->root_pos);
      pB = EulerSlerp(-1.0/3, endPosture->root_pos, an_);
    }
    else if (nextPosture == NULL)
    {
      an_ = EulerSlerp(1.0/2, EulerDouble(prevPosture->root_pos, startPosture->root_pos), endPosture->root_pos);
      pA = EulerSlerp(1.0/3, startPosture->root_pos, an_);

      pB = EulerSlerp(1.0/3, endPosture->root_pos, EulerDouble(prevPosture->root_pos, startPosture->root_pos));
    }
    else
    {
      an_ = EulerSlerp(1.0/2, EulerDouble(prevPosture->root_pos, startPosture->root_pos), endPosture->root_pos);
      pA = EulerSlerp(1.0/3, startPosture->root_pos, an_);

      an_ = EulerSlerp(1.0/2, EulerDouble(startPosture->root_pos, endPosture->root_pos), nextPosture->root_pos);
      pB = EulerSlerp(-1.0/3, endPosture->root_pos, an_);
    }

    pStart = startPosture->root_pos;
    pEnd = endPosture->root_pos;

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    // interpolate in between
    for(int frame=1; frame<=N; frame++)
    {
      Posture interpolatedPosture;
      Quaternion<double> q0, q1, q2, q3, interpolatedPostureQuaternion;
      Quaternion<double> temp, a, b;
      double t = 1.0 * frame / (N+1);

      // interpolate root position
      interpolatedPosture.root_pos = DeCasteljauEuler(t, pStart, pA, pB, pEnd);

      // interpolate bone rotations
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
      {
        if (startKeyframe - N - 1 > 0)
        {
          Euler2Quaternion(prevPosture->bone_rotation[bone].p, q0);
        }
        Euler2Quaternion(startPosture->bone_rotation[bone].p, q1);
        Euler2Quaternion(endPosture->bone_rotation[bone].p, q2);
        if (endKeyframe + N + 1 < inputLength)
        {
          Euler2Quaternion(nextPosture->bone_rotation[bone].p, q3);
        }

        if (prevPosture == NULL)
        {
          a = Slerp(1.0/3, q1, Double(q3, q2));

          temp = Slerp(0.5, Double(q1, q2), q3);
          b = Slerp(-1.0/3, q2, temp);
        }
        else if (nextPosture == NULL)
        {
          temp = Slerp(1.0/2, Double(q0, q1), q2);
          a = Slerp(1.0/3, q2, temp);

          b = Slerp(1.0/3, q2, Double(q0, q1));
        }
        else
        {
          temp = Slerp(1.0/2, Double(q0, q1), q2);
          a = Slerp(1.0/3, q1, temp);

          temp = Slerp(1.0/2, Double(q1, q2), q3);
          b = Slerp(-1.0/3, q2, temp);
        }

        interpolatedPostureQuaternion = DeCasteljauQuaternion(t, q1, a, b, q2);
        Quaternion2Euler(interpolatedPostureQuaternion, interpolatedPosture.bone_rotation[bone].p);
      }
      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;
  }

  for(int frame=startKeyframe+1; frame<inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}
