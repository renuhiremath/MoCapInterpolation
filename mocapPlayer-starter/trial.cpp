
void Interpolator::BezierInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

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

    vector pStart, pEnd;
    vector pA, pB, temp;
    vector p0, p1, p2, p3;

    // if (startKeyframe - N - 1 > 0)
    // {
    //   p0 = prevPosture->root_pos;
    // }
    // p1 = startPosture->root_pos;
    // p2 = endPosture->root_pos;
    // if (endKeyframe + N + 1 < inputLength)
    // {
    //   p3 = nextPosture->root_pos;
    // }
    //
    // //find the other control points for the spline
    // if (startKeyframe - N - 1 < 0) //p0 doesn't exist
    // {
    //   temp = p2 - p3 + p2;
    //   pA = (temp - p1) * (1.0/3) + p1;
    //
    //   temp = p2 - p1 + p2;
		// 	temp = (temp + p3) * 0.5;
		// 	temp = (temp - p2) * (1.0 / 3.0);
		// 	pB = p2 - temp;
    // }
    // else if (endKeyframe + N + 1 > inputLength) // p3 doesn't exist
    // {
    //   temp = p1 - p0 + p1;
    //   temp = (temp + p2) / 2.0;
    //   pA = (temp - p1) * (1.0 / 3.0) + p1;
    //
    //   temp = p1 - p0 + p1;
    //   pB = (temp - p2) * (1.0 / 3.0) + p2;
    // }
    // else
    // {
    //   temp = (((p1 - p0) + p1) + p2) / 2.0;
    //   pA = ((temp - p1) / 3.0) + p1;
    //
    //   temp = (((p2 - p1) + p2) + p3) / 2.0;
    //   temp = (temp - p2) * (1.0/3);
    //   pB = p2 - temp;
    // }
    //
    // pStart = p1;
    // pEnd = p2;

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    // interpolate in between
    for(int frame=1; frame<=N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N+1);
      vector a, b;

      // interpolate root position
      //interpolatedPosture.root_pos = DeCasteljauEuler(t, pStart, pA, pB, pEnd);
      interpolatedPosture.root_pos = startPosture->root_pos * (1-t) + endPosture->root_pos * t;

      // interpolate bone rotations
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
      {
        if (startKeyframe - N - 1 > 0)
        {
          p0 = prevPosture->bone_rotation[bone];
        }
        p1 = startPosture->bone_rotation[bone];
        p2 = endPosture->bone_rotation[bone];
        if (endKeyframe + N + 1 < inputLength)
        {
          p3 = nextPosture->bone_rotation[bone];
        }

        //find the other control points for the spline
        if (startKeyframe - N - 1 < 0) //p0 doesn't exist
        {
          temp = p2 - p3 + p2;
          a = (temp - p1) * (1.0/3) + p1;

          temp = p2 - p1 + p2;
    			temp = (temp + p3) * 0.5;
    			temp = (temp - p2) * (1.0 / 3.0) + p2;
    			b = p2 - temp + p2;
        }
        else if (endKeyframe + N + 1 > inputLength) // p3 doesn't exist
        {
          temp = p1 - p0 + p1;
          temp = (temp + p2) * 0.5;
          a = (temp - p1) * (1.0 / 3.0) + p1;

          temp = p1 - p0 + p1;
          b = (temp - p2) * (1.0 / 3.0) + p2;
        }
        else
        {
          temp = (((p1 - p0) + p1) + p2) * 0.5;
          a = (temp - p1) * (1.0/3) + p1;

          temp = (((p2 - p1) + p2) +p3) * 0.5;
          temp = (temp - p2) * (1.0/3) + p2;
          b = p2 - temp + p2;
        }

        interpolatedPosture.bone_rotation[bone] = DeCasteljauEuler(t, p1, a, b, p2);
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

  int startKeyframe = 0;
  while (startKeyframe + N + 1 < inputLength)
  {
    int endKeyframe = startKeyframe + N + 1;
    Posture * prevPosture;
    Posture * nextPosture;

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

    Quaternion<double> q0, q1, q2, q3, pStart, pEnd;
    Quaternion<double> pA, pB, temp;
    Quaternion<double> p0, p1, p2, p3;

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    // interpolate in between
    for(int frame=1; frame<=N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N+1);
      Quaternion<double> a, b, interpolatedPostureQuaternion;

      // interpolate root position
      interpolatedPosture.root_pos = startPosture->root_pos * (1-t) + endPosture->root_pos * t;

      // interpolate bone rotations
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
      {
        if (startKeyframe - N - 1 > 0)
        {
          Euler2Quaternion(prevPosture->bone_rotation[bone].p, p0);
        }
        Euler2Quaternion(startPosture->bone_rotation[bone].p, p1);
        Euler2Quaternion(endPosture->bone_rotation[bone].p, p2);
        if (endKeyframe + N + 1 < inputLength)
        {
          Euler2Quaternion(nextPosture->bone_rotation[bone].p, p3);
        }

        //find the other control points for the spline
        if (startKeyframe - N - 1 < 0) //p0 doesn't exist
        {
          temp = p2 - p3 + p2;
          a = (temp - p1) * (1.0/3) + p1;

          temp = p2 - p1 + p2;
    			temp = (temp + p3) * 0.5;
    			temp = (temp - p2) * (1.0 / 3.0) + p2;
    			b = p2 - temp + p2;
        }
        else if (endKeyframe + N + 1 > inputLength) // p3 doesn't exist
        {
          temp = p1 - p0 + p1;
          temp = (temp + p2) * 0.5;
          a = (temp - p1) * (1.0 / 3.0) + p1;

          temp = p1 - p0 + p1;
          b = (temp - p2) * (1.0 / 3.0) + p2;
        }
        else
        {
          temp = (((p1 - p0) + p1) + p2) * 0.5;
          a = (temp - p1) * (1.0/3) + p1;

          temp = (((p2 - p1) + p2) +p3) * 0.5;
          temp = (temp - p2) * (1.0/3) + p2;
          b = p2 - temp + p2;
        }
        interpolatedPostureQuaternion = DeCasteljauQuaternion(t, p1, a, b, p2);
        Quaternion2Euler(interpolatedPostureQuaternion, interpolatedPosture.bone_rotation[bone].p);
      }

      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;
  }

  for(int frame=startKeyframe+1; frame<inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}
