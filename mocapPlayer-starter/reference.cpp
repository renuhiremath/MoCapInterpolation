{
/*
  interpolator.h

  Create interpolated motion.

  Revision 1 - Alla and Kiran, Jan 18, 2002
  Revision 2 - Jernej Barbic, Yili Zhao, Feb 2012
*/


#ifndef _INTERPOLATOR_H
#define _INTERPOLATOR_H


enum InterpolationType
{
  LINEAR = 0, BEZIER = 1
};

enum AngleRepresentation
{
  EULER = 0, QUATERNION = 1
};

class Interpolator
{
public:
  //constructor, destructor
  Interpolator();
  ~Interpolator();

  //Set interpolation type
  void SetInterpolationType(InterpolationType interpolationType) {m_InterpolationType = interpolationType;};
  //Set angle representation for interpolation
  void SetAngleRepresentation(AngleRepresentation angleRepresentation) {m_AngleRepresentation = angleRepresentation;};

  //Create interpolated motion and store it into pOutputMotion (which will also be allocated)
  void Interpolate(Motion * pInputMotion, Motion ** pOutputMotion, int N);

protected:
  InterpolationType m_InterpolationType; //Interpolation type (Linear, Bezier)
  AngleRepresentation m_AngleRepresentation; //Angle representation (Euler, Quaternion)

  // conversion routines
  // angles are given in degrees; assume XYZ Euler angle order
  void Rotation2Euler(double R[9], double angles[3]);
  void Euler2Rotation(double angles[3], double R[9]);
  void Euler2Quaternion(double angles[3], Quaternion<double> & q);
  void Quaternion2Euler(Quaternion<double> & q, double angles[3]);

  // quaternion interpolation
  Quaternion<double> Slerp(double t, Quaternion<double> & qStart, Quaternion<double> & qEnd);
  Quaternion<double> Double(Quaternion<double> p, Quaternion<double> q);

  // interpolation routines
  void LinearInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N);
  void BezierInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N);
  void LinearInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N);
  void BezierInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N);

  // Bezier spline evaluation
  vector DeCasteljauEuler(double t, vector p0, vector p1, vector p2, vector p3); // evaluate Bezier spline at t, using DeCasteljau construction, vector version
  Quaternion<double> DeCasteljauQuaternion(double t, Quaternion<double> p0, Quaternion<double> p1, Quaternion<double> p2, Quaternion<double> p3); // evaluate Bezier spline at t, using DeCasteljau construction, Quaternion version

};

#endif
}
{
  /*
  motion.h

  Motion class

  1. read an AMC file and store it in a sequence of state vector
  2. write an AMC file
  3. export to a mrdplot format for plotting the trajectories

  You can add more motion data processing functions in this class.

  Revision 1 - Steve Lin, Jan. 14, 2002
  Revision 2 - Alla and Kiran, Jan 18, 2002
  Revision 3 - Jernej Barbic and Yili Zhao, Feb 2012
  */

  #ifndef _MOTION_H_
  #define _MOTION_H_

  #include "vector.h"
  #include "types.h"
  #include "posture.h"
  #include "skeleton.h"

  class Motion
  {
    //function members
  public:

    // parse AMC file (default scale=0.06)
    Motion(char *amc_filename, double scale, Skeleton * pSkeleton);

    //Use to create default motion with specified number of frames
    Motion(int numFrames, Skeleton * pSkeleton);

    ~Motion();

    // scale is a parameter to adjust the translationalal scaling
    // the value of scale should be consistent with the scale parameter used in Skeleton()
    // forceAllJointsBe3DOF should be set to 0; use 1 to signal that the file contains three Euler
    // angles for all the joints, even those that are 1-dimensional or 2-dimensional (advanced usage)
    int writeAMCfile(char* filename, double scale, int forceAllJointsBe3DOF=0);

    //Set all postures to default posture
    //Root position at (0,0,0), orientation of each bone to (0,0,0)
    void SetPosturesToDefault();

    //Set the entire posture at specified frame (posture = root position and all bone rotations)
    void SetPosture(int frameIndex, Posture InPosture);

    //Set root position at specified frame
    void SetRootPos(int frameIndex, vector vPos);
    //Set specified bone rotation at specified frame
    void SetBoneRotation(int frameIndex, int boneIndex, vector vRot);

    int GetNumFrames() { return m_NumFrames; }
    Posture * GetPosture(int frameIndex);

    Skeleton * GetSkeleton() { return pSkeleton; }

  protected:
    int m_NumFrames; //number of frames in the motion
    Skeleton * pSkeleton;
    //Root position and all bone rotation angles for each frame (as read from AMC file)
    Posture * m_pPostures;

    // The default value is 0.06
    int readAMCfile(char* name, double scale);
  };

  #endif


}
{
  #ifndef _QUATERNION_H_
  #define _QUATERNION_H_

  /*
  Quaternion C++ class.
  This class implements quaternions and the commonly used
  algebraic operations on quaternions. The class is templated:
  you can use either float or double precision.
  Supports using quaternions to represent/manipulate rotations.

  q = s + x * i + y * j + z * k

  If using quaternions to represent rotations, q must be a unit quaternion.
    Then the following relates q to the corresponding rotation:
      s = cos(angle/2)
      (x,y,z) = sin(angle/2) * axis_of_rotation
      (axis_of_rotation is unit length)

  See also example.cpp .

  Version: 1.0
  */

  #include <stdio.h>
  #include <math.h>

  // forward declarations for the friend template
  template <typename real> class Quaternion;
  template <typename real> Quaternion<real> operator* (real alpha, Quaternion<real> q2);

  template <typename real>
  class Quaternion
  {
  public:

    inline Quaternion(); // q = 0
    inline Quaternion(real s, real x, real y, real z); // q = s + x * i + y * j + z * k
    inline Quaternion(real s); // q = s + 0 * i + 0 * j + 0 * k;

    // Makes the unit quaternion corresponding to a rotation around axis 'unitAxis' of angle 'angle'.
    // Angle is in radians; unitAxis must be a unit 3D vector.
    inline Quaternion(real angle, real unitAxis[3]);

    inline void Set(real s, real x, real y, real z); // sets quaternion to the new value

    inline real Gets() const;
    inline real Getx() const;
    inline real Gety() const;
    inline real Getz() const;

    inline Quaternion operator+ (const Quaternion q2) const; // q3 = q1+q2
    inline Quaternion operator- (const Quaternion q2) const; // q3 = q1-q2
    inline Quaternion operator* (const Quaternion q2) const; // q3 = q1 * q2
    inline Quaternion operator/ (const Quaternion q2) const; // q3 = q1 / q2

    // Multiply quaternion with a scalar; e.g. q1 = alpha * q2;
    friend Quaternion<real> operator* (real alpha, const Quaternion<real> q2)
    {
      return Quaternion<real>(alpha * q2.s, alpha * q2.x, alpha * q2.y, alpha * q2.z);
    }

    inline Quaternion conj(); // q2 = q1.conj()

    inline Quaternion & operator= (const Quaternion rhs); // q2 = q1;
    inline Quaternion & operator= (real s); // sets quaternion equal to the scalar quaternion s
    inline int operator== (const Quaternion rhs) const; // q2 == q1
    inline int operator!= (const Quaternion rhs) const; // q2 != q1

    void Normalize(); // q.Normalize() scales q such that it is unit size

    inline void MoveToRightHalfSphere(); //  if scalar part (that is, 's') is negative, this will multiply the quaternion by -1

    inline real Norm2() const; // returns the squared norm of the quaternion, i.e. s*s + x*x + y*y + z*z
    inline real Norm() const { return sqrt(Norm2()); }

    // Transforms the quaternion to the corresponding rotation matrix.
    // Quaternion is assumed to be a unit quaternion.
    // R is a 3x3 orthogonal matrix and will be returned in row-major order.
    inline void Quaternion2Matrix(real * R) const;

    // Transforms the given matrix (assumed orthogonal) into one of the two corresponding quaternions.
    // Matrix is assumed to be in row-major order.
    // There are two quaternions corresponding to a rotation (and they have opposite signs). You can't directly control which one you get, but you can force the real part to be non-negative by a subsequent call to MoveToRightHalfSphere() .
    // This implementation follows David Baraff's SIGGRAPH course notes:
    // http://www.cs.cmu.edu/~baraff/pbm/pbm.html
    static Quaternion Matrix2Quaternion(real * R);

    // Returns the angle of rotation (in radians), and the unit rotation axis corresponding to the quaternion.
    // Assumes a unit quaternion (use Normalize() to remove any noise due to floating point errors).
    // If s >= 0, the angle will be on the interval [0,pi] .
    // If s < 0, the angle will be on the interval (pi,2pi]. To get a representation where the angle is on [0, pi], you can manually flip the sign of the returned unitAxis, and use the angle of 2pi-angle. Alternatively, you can use MoveToRightHalfSphere before calling GetRotation to ensure that you are always in the s >= 0 case.
    inline void GetRotation(real * angle, real unitAxis[3]) const;

    // Returns (x,y,z) = sin(theta/2) * axis, where
    //   theta is the angle of rotation, theta is on [-pi,pi), and axis is the unit axis of rotation.
    // Assumes a unit quaternion.
    // Note: this routine is a bit exotic; I expect it to be not so widely used.
    inline void GetSinExponential(real * x, real * y, real * z) const;

    // Prints the quaternion to stdout.
    inline void Print() const;

  protected:
    real s,x,y,z;

  };

  template <typename real>
  inline Quaternion<real>::Quaternion()
  {
    s = x = y = z = 0;
  }

  template <typename real>
  inline Quaternion<real>::Quaternion(real s_)
  {
    s = s_;
    x = y = z = 0;
  }

  template <typename real>
  inline Quaternion<real>::Quaternion(real angle, real unitAxis[3])
  {
    s = cos(angle/2.0);
    real sin2 = sin(angle/2.0);
    x = sin2 * unitAxis[0];
    y = sin2 * unitAxis[1];
    z = sin2 * unitAxis[2];
  }

  template <typename real>
  inline Quaternion<real>::Quaternion(real s_, real x_, real y_, real z_)
  {
    s = s_;
    x = x_;
    y = y_;
    z = z_;
  }

  template <typename real>
  inline void Quaternion<real>::Set(real s_g, real x_g, real y_g, real z_g) // sets quaternion to the new value
  {
    s = s_g;
    x = x_g;
    y = y_g;
    z = z_g;
  }

  template <typename real>
  inline real Quaternion<real>::Gets() const { return s; }

  template <typename real>
  inline real Quaternion<real>::Getx() const { return x; }

  template <typename real>
  inline real Quaternion<real>::Gety() const { return y; }

  template <typename real>
  inline real Quaternion<real>::Getz() const { return z; }

  template <typename real>
  inline Quaternion<real> & Quaternion<real>::operator= (const Quaternion<real> rhs)
  {
    s = rhs.s;
    x = rhs.x;
    y = rhs.y;
    z = rhs.z;

    return *this;
  }

  template <typename real>
  inline Quaternion<real> & Quaternion<real>::operator= (real s_g)
  {
    s = s_g;
    x = 0;
    y = 0;
    z = 0;

    return *this;
  }

  template <typename real>
  inline int Quaternion<real>::operator== (const Quaternion<real> rhs) const
  {
    return ((s == rhs.s) && (x == rhs.x) &&
            (y == rhs.y) && (z == rhs.z));
  }

  template <typename real>
  inline int Quaternion<real>::operator!= (const Quaternion<real> rhs) const
  {
    return ((s != rhs.s) || (x != rhs.x) ||
            (y != rhs.y) || (z != rhs.z));
  }

  template <typename real>
  inline void Quaternion<real>::Normalize()
  {
    real invNorm;
    invNorm = (real)1.0 / (real)sqrt(Norm2());

    s *= invNorm;
    x *= invNorm;
    y *= invNorm;
    z *= invNorm;

  }

  template <typename real>
  inline real Quaternion<real>::Norm2() const
  {
    return (s*s + x*x + y*y + z*z);
  }

  template <typename real>
  inline Quaternion<real> Quaternion<real>::operator+ (const Quaternion<real> q2) const
  {
    Quaternion<real> w(s + q2.s, x + q2.x, y + q2.y, z + q2.z);

    return w;
  }

  template <typename real>
  inline Quaternion<real> Quaternion<real>::operator- (const Quaternion<real> q2) const
  {
    Quaternion<real> w(s - q2.s, x - q2.x, y - q2.y, z - q2.z);
    return w;
  }

  template <typename real>
  inline Quaternion<real> Quaternion<real>::operator* (const Quaternion<real> q2) const
  {
    Quaternion<real> w(
          s * q2.s - x * q2.x - y    * q2.y - z * q2.z,
          s * q2.x + q2.s * x + y    * q2.z - q2.y * z,
          s * q2.y + q2.s * y + q2.x * z    - x    * q2.z,
          s * q2.z + q2.s * z + x    * q2.y - q2.x * y);

    return w;
  }

  template <typename real>
  inline Quaternion<real> Quaternion<real>::operator/ (const Quaternion<real> q2) const
  {
    // compute invQ2 = q2^{-1}
    Quaternion<real> invQ2;
    real invNorm2 = 1.0 / q2.Norm2();
    invQ2.s = q2.s * invNorm2;
    invQ2.x = -q2.x * invNorm2;
    invQ2.y = -q2.y * invNorm2;
    invQ2.z = -q2.z * invNorm2;

    // result = *this * invQ2
    return (*this * invQ2);

  }

  template <typename real>
  inline Quaternion<real> Quaternion<real>::conj()
  {
    Quaternion<real> w(s,-x,-y,-z);
    return w;
  }

  // Transforms the quaternion to the corresponding rotation matrix.
  // Quaternion is assumed to be a unit quaternion.
  // R is a 3x3 orthogonal matrix and will be returned in row-major order.
  template <typename real>
  inline void Quaternion<real>::Quaternion2Matrix(real * R) const
  {
    R[0] = 1 - 2*y*y - 2*z*z; R[1] = 2*x*y - 2*s*z;     R[2] = 2*x*z + 2*s*y;
    R[3] = 2*x*y + 2*s*z;     R[4] = 1 - 2*x*x - 2*z*z; R[5] = 2*y*z - 2*s*x;
    R[6] = 2*x*z - 2*s*y;     R[7] = 2*y*z + 2*s*x;     R[8] = 1 - 2*x*x - 2*y*y;
  }

  // Returns (x,y,z) = sin(theta/2) * axis, where
  //   theta is the angle of rotation, theta\in\{-pi,pi\}, and
  //   axis is the unit axis of rotation.
  template <typename real>
  inline void Quaternion<real>::GetSinExponential(real * sex, real * sey, real * sez)const
  {
    if (s<0)
    {
      *sex = -x;
      *sey = -y;
      *sez = -z;
    }
    else
    {
      *sex = x;
      *sey = y;
      *sez = z;
    }
  }

  template <typename real>
  inline void Quaternion<real>::GetRotation(real * angle, real unitAxis[3]) const
  {
    if ((s >= ((real)1)) || (s <= (real)(-1)))
    {
      // identity; this check is necessary to avoid problems with acos if s is 1 + eps
      *angle = 0;
      unitAxis[0] = 1;
      unitAxis[0] = 0;
      unitAxis[0] = 0;
      return;
    }

    *angle = 2.0 * acos(s);
    real sin2 = x*x + y*y + z*z; //sin^2(*angle / 2.0)

    if (sin2 == 0)
    {
      // identity rotation; angle is zero, any axis is equally good
      unitAxis[0] = 1;
      unitAxis[0] = 0;
      unitAxis[0] = 0;
    }
    else
    {
      real inv = 1.0 / sqrt(sin2); // note: *angle / 2.0 is on [0,pi], so sin(*angle / 2.0) >= 0, and therefore the sign of sqrt can be safely taken positive
      unitAxis[0] = x * inv;
      unitAxis[1] = y * inv;
      unitAxis[2] = z * inv;
    }
  }

  template <typename real>
  inline void Quaternion<real>::MoveToRightHalfSphere()
  {
    if (s<0)
    {
      s *= -1;
      x *= -1;
      y *= -1;
      z *= -1;
    }
  }

  template <typename real>
  inline void Quaternion<real>::Print() const
  {
    printf("%f + %fi + %fj + %fk\n",s,x,y,z);
  }

  #endif

}
/*

Revision 1 - Steve Lin, Jan. 14, 2002
Revision 2 - Alla and Kiran, Jan 18, 2002
Revision 3 - Jernej Barbic and Yili Zhao, Feb, 2012

*/
#ifndef _TYPES_H
#define _TYPES_H

// Use this parameter to adjust the size of the skeleton. The default value is 0.06.
// This creates a human skeleton of 1.7 m in height (approximately)
//#define MOCAP_SCALE 0.06
#define MOCAP_SCALE 0.06
//static const int	NUM_BONES_IN_ASF_FILE	= 31;
#define MAX_BONES_IN_ASF_FILE 256
#define MAX_CHAR 1024
#define MAX_SKELS 16

#define PM_MAX_FRAMES 60000

#ifndef M_PI
#define M_PI 3.14159265
#endif

enum ErrorType
{
  NO_ERROR_SET = 0, BAD_OFFSET_FILE, NOT_SUPPORTED_INTERP_TYPE, BAD_INPUT_FILE
};


#endif
/*
display.h

Display the skeleton, ground plane and other objects.

Revision 1 - Steve Lin, Jan. 14, 2002
Revision 2 - Alla and Kiran, Jan 18, 2002
Revision 3 - Jernej Barbic and Yili Zhao, Feb, 2012
*/

#ifndef _DISPLAY_SKELETON_H_
#define _DISPLAY_SKELETON_H_

#include <FL/glu.h>
#include "skeleton.h"
#include "motion.h"

class DisplaySkeleton
{

  //member functions
public:
  enum RenderMode
  {
    BONES_ONLY, BONES_AND_LOCAL_FRAMES
  };
  enum JointColor
  {
    GREEN, RED, BLUE, NUMBER_JOINT_COLORS
  };

  DisplaySkeleton();
  ~DisplaySkeleton();

  //set skeleton for display
  void LoadSkeleton(Skeleton * pSkeleton);
  //set motion for display
  void LoadMotion(Motion * pMotion);

  //display the scene (skeleton, ground plane ....)
  void Render(RenderMode renderMode);
  void RenderShadow(double ground[4], double light[4]);

  void SetDisplayedSpotJoint(int jointID) {m_SpotJoint = jointID;}
  int GetDisplayedSpotJoint(void) {return m_SpotJoint;}
  int GetNumSkeletons(void) {return numSkeletons;}
  Skeleton * GetSkeleton(int skeletonIndex);
  Motion * GetSkeletonMotion(int skeletonIndex);

  void Reset(void);

protected:
  RenderMode renderMode;
  // Draw a particular bone
  void DrawBone(Bone *ptr, int skelNum);
  // Draw the skeleton hierarchy
  void Traverse(Bone *ptr, int skelNum);
  // Model matrix for the shadow
  void SetShadowingModelviewMatrix(double ground[4], double light[4]);
  void DrawSpotJointAxis(void);
  void SetDisplayList(int skeletonID, Bone *bone, GLuint *pBoneList);

  int m_SpotJoint;		//joint whose local coordinate framework is drawn
  int numSkeletons;
  Skeleton *m_pSkeleton[MAX_SKELS];		//pointer to current skeleton
  Motion *m_pMotion[MAX_SKELS];		//pointer to current motion
  GLuint m_BoneList[MAX_SKELS];		//display list with bones

  static float jointColors[NUMBER_JOINT_COLORS][3];
};

#endif
// generated by Fast Light User Interface Designer (fluid) version 1.0009
/*
Revision 1 - Steve Lin, Jan. 14, 2002
Revision 2 - Alla and Kiran, Jan 18, 2002
Revision 3 - Jernej Barbic and Yili Zhao, Feb, 2012
*/

#ifndef interface_h
#define interface_h

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Value_Input.H>
#include <FL/Fl_Light_Button.H>
#include <FL/Fl_Value_Slider.H>
#include "mocapPlayer.h"

extern Fl_Window * main_window;

extern Fl_Button * play_button;
extern Fl_Button * pause_button;
extern Fl_Button * rewind_button;
extern Fl_Button * repeat_button;
extern Fl_Button * plusOne_button;
extern Fl_Button * minusOne_button;
extern void play_callback(Fl_Button *, void *);

extern Fl_Button * loadSkeleton_button;
extern Fl_Button * loadMotion_button;
extern void load_callback(Fl_Button *, void *);

extern Fl_Button * reloadMotion_button;
extern void reload_callback(Fl_Button *, void *);

extern Fl_Button * resetScene_button;
extern void resetScene_callback(Fl_Button *, void *);

extern Fl_Button * screenShot_button;
extern void saveScreenToFile_callback(Fl_Button *, void *);

extern Fl_Light_Button * record_button;
extern void record_callback(Fl_Light_Button *, void *);

extern Fl_Value_Slider * frame_slider;
extern void fslider_callback(Fl_Value_Slider * , long);

extern Fl_Light_Button * groundPlane_button;
extern void renderGroundPlane_callback(Fl_Light_Button *, long);


extern Fl_Light_Button * worldAxes_button;
extern void renderWorldAxes_callback(Fl_Light_Button *, long);

extern Fl_Light_Button * fog_button;
extern void useFog_callback(Fl_Light_Button *, long);

extern Fl_Button * aboutPlayer_button;
extern void aboutPlayer_callback(Fl_Button * button, void *);

extern Player_Gl_Window *glwindow;

extern Fl_Value_Input * speedUp;
extern void playSpeed_callback(Fl_Value_Input *, void *);

extern Fl_Value_Input * joint_idx;
extern void spotJoint_callback(Fl_Value_Input *, void *);

extern Fl_Value_Input * sub_input;
extern void skeletonID_callback(Fl_Value_Input * , void *);

extern Fl_Value_Input * tx_input;
extern void tx_callback(Fl_Value_Input *, void *);

extern Fl_Value_Input * ty_input;
extern void ty_callback(Fl_Value_Input *, void *);

extern Fl_Value_Input * tz_input;
extern void tz_callback(Fl_Value_Input *, void *);

extern Fl_Value_Input * rx_input;
extern void rx_callback(Fl_Value_Input *, void *);

extern Fl_Value_Input * ry_input;
extern void ry_callback(Fl_Value_Input*, void *);

extern Fl_Value_Input * rz_input;
extern void rz_callback(Fl_Value_Input *, void *);

Fl_Window * make_window();

#endif
/*
Revision 1 - Steve Lin, Jan. 14, 2002
Revision 2 - Alla and Kiran, Jan 18, 2002
Revision 3 - Jernej Barbic and Yili Zhao, Feb, 2012
*/
#ifndef _POSTURE_H
#define _POSTURE_H

#include "vector.h"
#include "types.h"

//Root position and all bone rotation angles (including root)
struct Posture
{
public:
  //Root position (x, y, z)
  vector root_pos;

  //Euler angles (thetax, thetay, thetaz) of all bones in their local coordinate system.
  //If a particular bone does not have a certain degree of freedom,
  //the corresponding rotation is set to 0.
  //The order of the bones in the array corresponds to their ids in .ASf file: root, lhipjoint, lfemur, ...
  vector bone_rotation[MAX_BONES_IN_ASF_FILE];

  // bones that are translated relative to parents (resulting in gaps) (rarely used)
  vector bone_translation[MAX_BONES_IN_ASF_FILE];

  // bones that change length during the motion (rarely used)
  vector bone_length[MAX_BONES_IN_ASF_FILE];
};

#endif
#ifndef _PIC_H
#define _PIC_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

#ifndef ALLOC
#define ALLOC(ptr, type, n) { \
    ptr = (type *)malloc((n)*sizeof(type)); \
    if (!ptr) { \
	fprintf(stderr, "pic_alloc: Can't allocate %d bytes of memory, aborting\n", n); \
	exit(1); \
    } \
}
#endif

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

typedef unsigned char Pixel1;	/* 1-byte pixel */
typedef struct {Pixel1 r, g, b;}	Pixel1_rgb;
typedef short                           Pixel2;
typedef struct {Pixel2 r, g, b;}        Pixel2_rgb;
typedef Pixel2_rgb Rgbcolor;

typedef struct {		/* PICTURE */
    int nx, ny;			/* width & height, in pixels */
    int bpp;			/* bytes per pixel = 1, 3, or 4 */
    Pixel1 *pix;		/* array of pixels */
				/* data is in row-major order,
				    i.e. it has the same memory layout as:
				    if 1 byte per pixel,  then array[ny][nx]
				    if 3 bytes per pixel, then array[ny][nx][3]
				    if 4 bytes per pixel, then array[ny][nx][4] */
} Pic;

/*
 * use the following macro to access the red, green, and blue channels of
 * a given pixel, for an image stored in a Pic structure.
 */
#define PIC_PIXEL(pic, x, y, chan) \
    (pic)->pix[((y)*(pic)->nx+(x))*(pic)->bpp+(chan)]
    /* returns channel chan of pixel (x,y) of picture pic */
    /* usually chan=0 for red, 1 for green, 2 for blue */

typedef enum {PIC_TIFF_FILE, PIC_PPM_FILE, PIC_JPEG_FILE, PIC_UNKNOWN_FILE} Pic_file_format; // only PPM is supported

/*----------------------Allocation routines--------------------------*/
extern Pic *pic_alloc(int nx, int ny, int bytes_per_pixel, Pic *opic);
extern void pic_free(Pic *p);

/*------------------------- I/O routines ----------------------------*/
/*
extern int tiff_get_size(char *file, int *nx, int *ny);
extern Pic *tiff_read(char *file, Pic *opic);
extern int tiff_write(char *file, Pic *pic);

extern int jpeg_get_size(char *file, int *nx, int *ny);
extern Pic *jpeg_read(char *file, Pic *opic);
extern int jpeg_write(char *file, Pic *pic);
*/

extern int ppm_get_size(char *file, int *nx, int *ny);
extern Pic *ppm_read(char *file, Pic *opic);
extern int ppm_write(char *file, Pic *pic);

/*
extern int pic_get_size(char *file, int *nx, int *ny);
extern Pic *pic_read(char *file, Pic *opic);
extern int pic_write(char *file, Pic *pic, Pic_file_format format);
extern Pic_file_format pic_file_type(char *file);
extern Pic_file_format pic_filename_type(char *file);
*/

#ifdef __cplusplus
}
#endif
#endif
/*
player.h

All the user interface functions .

Revision 1 - Steve Lin, Jan. 14, 2002
Revision 2 - Alla and Kiran, Jan 18, 2002
Revision 3 - Jernej Barbic and Yili Zhao, Feb, 2012
*/

#ifndef _PLAYER_H
#define _PLAYER_H

#include <FL/Fl_Gl_Window.H>

class Player_Gl_Window : public Fl_Gl_Window
{
public:
  inline Player_Gl_Window(int x, int y, int w, int h, const char *l=0) :
  Fl_Gl_Window(x, y, w, h, l) {};

  /* This is an overloading of a Fl_Gl_Window call.  It is called
  whenever a window needs refreshing. */
  void draw();

  /* This is an overloading of a Window call.  It is
  called whenever a event happens inside the space
  taken up by the Anim_Gl_Window. */
  int handle(int event);
};


typedef struct _MouseT {
  int button;
  int x;
  int y;
} MouseT;


typedef struct _CameraT {
  double zoom;
  double tw;
  double el;
  double az;
  double tx;
  double ty;
  double tz;
  double atx;
  double aty;
  double atz;
} CameraT;


void GraphicsInit();
void display();

#endif
/*
skeleton.h

Definition of the skeleton.

Revision 1 - Steve Lin, Jan. 14, 2002
Revision 2 - Alla and Kiran, Jan 18, 2002
Revision 3 - Jernej Barbic and Yili Zhao, Feb, 2012

*/

#ifndef _SKELETON_H
#define _SKELETON_H

#include "posture.h"

// this structure defines the property of each bone segment, including its connection to other bones,
// DOF (degrees of freedom), relative orientation and distance to the outboard bone
struct Bone
{
  struct Bone *sibling;	// Pointer to the sibling (branch bone) in the hierarchy tree
  struct Bone *child; // Pointer to the child (outboard bone) in the hierarchy tree

  int idx; // Bone index

  double dir[3]; // Unit vector describes the direction from local origin to
  // the origin of the child bone
  // Notice: stored in local coordinate system of the bone

  double length; // Bone length

  double axis_x, axis_y, axis_z;// orientation of each bone's local coordinate
  //system as specified in ASF file (axis field)

  double aspx, aspy; // aspect ratio of bone shape

  int dof; // number of bone's degrees of freedom
  int dofrx, dofry, dofrz; // rotational degree of freedom mask in x, y, z axis
  int doftx, dofty, doftz; // translational degree of freedom mask in x, y, z axis
  int doftl;
  // dofrx=1 if this bone has x rotational degree of freedom, otherwise dofrx=0.

  // bone names
  char name[256];
  // rotation matrix from the local coordinate of this bone to the local coordinate system of it's parent
  double rot_parent_current[4][4];

  //Rotation angles for this bone at a particular time frame (as read from AMC file) in local coordinate system,
  //they are set in the setPosture function before display function is called
  double rx, ry, rz;
  double tx,ty,tz;
  double tl;
  int dofo[8];
};


class Skeleton
{
public:
  // The scale parameter adjusts the size of the skeleton. The default value is 0.06 (MOCAP_SCALE).
  // This creates a human skeleton of 1.7 m in height (approximately)
  Skeleton(char *asf_filename, double scale);
  ~Skeleton();

  //Get root node's address; for accessing bone data
  Bone* getRoot();
  static int getRootIndex() { return 0; }

  //Set the skeleton's pose based on the given posture
  void setPosture(Posture posture);

  //Initial posture Root at (0,0,0)
  //All bone rotations are set to 0
  void setBasePosture();

  // marks previously unavailable rotational DOFs as available, and sets them to 0
  void enableAllRotationalDOFs();

  int name2idx(char *);
  char * idx2name(int);
  void GetRootPosGlobal(double rootPosGlobal[3]);
  void GetTranslation(double translation[3]);
  void GetRotationAngle(double rotationAngle[3]);
  void SetTranslationX(double tx_){tx = tx_;}
  void SetTranslationY(double ty_){ty = ty_;}
  void SetTranslationZ(double tz_){tz = tz_;}
  void SetRotationAngleX(double rx_){rx = rx_;}
  void SetRotationAngleY(double ry_){ry = ry_;}
  void SetRotationAngleZ(double rz_){rz = rz_;}

  int numBonesInSkel(Bone bone);
  int movBonesInSkel(Bone bone);

protected:

  //parse the skeleton (.ASF) file
  int readASFfile(char* asf_filename, double scale);

  //This recursive function traverses skeleton hierarchy
  //and returns a pointer to the bone with index - bIndex
  //ptr should be a pointer to the root node
  //when this function first called
  Bone* getBone(Bone *ptr, int bIndex);

  //This function sets sibling or child for parent bone
  //If parent bone does not have a child,
  //then pChild is set as parent's child
  //else pChild is set as a sibling of parents already existing child
  int setChildrenAndSibling(int parent, Bone *pChild);

  //Rotate all bone's direction vector (dir) from global to local coordinate system
  void RotateBoneDirToLocalCoordSystem();

  void set_bone_shape(Bone *bone);
  void compute_rotation_parent_child(Bone *parent, Bone *child);
  void ComputeRotationToParentCoordSystem(Bone *bone);

  // root position in world coordinate system
  double m_RootPos[3];
  double tx,ty,tz;
  double rx,ry,rz;

  int NUM_BONES_IN_ASF_FILE;
  int MOV_BONES_IN_ASF_FILE;

  Bone *m_pRootBone;  // Pointer to the root bone, m_RootBone = &bone[0]
  Bone  m_pBoneList[MAX_BONES_IN_ASF_FILE];   // Array with all skeleton bones

  void removeCR(char * str); // removes CR at the end of line
};

#endif
/*
transform.h

Revision 1 - Steve Lin, Jan. 14, 2002
Revision 2 - Alla and Kiran, Jan 18, 2002
Revision 3 - Jernej Barbic and Yili Zhao, Feb, 2012

*/

#ifndef _TRANSFORM_H
#define _TRANSFORM_H

class Matrix
{

};

void matrix_transpose(double a[4][4], double b[4][4]);
void matrix_print(char *str, double a[4][4]);
void matrix_transform_affine(double m[4][4], double x, double y, double z, double pt[3]);
void matrix_mult(double a[][4], double b[][4], double c[][4]);
void matrix_mult_3x3(double a[][3], double b[][3], double c[][3]);

void v3_cross(double a[3], double b[3], double c[3]);
double v3_mag(double a[3]);
double v3_dot(double a[3], double b[3]);

//Rotate vector v around axis X by angle a, around axis Y by angle b and around axis Z by angle c
void vector_rotationXYZ(double *v, double a, double b, double c);

//Create Rotation matrix, that rotates around axis X by angle a
void rotationX(double r[][4], double a);
//Create Rotation matrix, that rotates around axis Y by angle a
void rotationY(double r[][4], double a);
//Create Rotation matrix, that rotates around axis Z by angle a
void rotationZ(double r[][4], double a);

//Return the angle between vectors v1 and v2 around the given axis
double GetAngle(double* v1, double* v2, double* axis);

#endif
/*

* Copyright (c) 2008, Carnegie Mellon University
* All rights reserved.
*
* Code author: Jernej Barbic
* Research: Jernej Barbic, Doug L. James
* Funding: NSF, Link Foundation
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of Carnegie Mellon University, nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


=== A counter to measure code execution time. ===

Designed for real-time system (e.g., real-time computer animation, haptics),
but useful in general. You can time arbitrary segments of your code.
Same interface under Windows, Linux and Mac OS X.

Under Linux/MAC OS X, accuracy is that of gettimeofday function,
which gives time in seconds and microseconds.
In practice, it has been accurate down to microsecond range.

Under Windows, the counter uses the QueryPerformanceCounter Windows API call.
Again, accuracy has been in the microsecond range in practice.

Usage:
Call StartCounter() before your code block.
Call StopCounter() after your code block.
Read the elapsed time using GetElapsedTime().

Version: 1.0

*/

#ifndef _PERFORMANCECOUNTER_H_
#define _PERFORMANCECOUNTER_H_

/**************** LINUX/MAC OS X COUNTER *******************/

#if (defined __unix__) || (defined __APPLE__)

#include "stdlib.h"
#include "sys/time.h"

class PerformanceCounter
{
public:

  PerformanceCounter() { StartCounter(); }

  void StartCounter(); // call this before your code block
  void StopCounter(); // call this after your code block

  // read elapsed time (units are seconds, accuracy is up to microseconds)
  double GetElapsedTime();

protected:
  long startCountSec,stopCountSec,startCountMicroSec,stopCountMicroSec;
};

inline void PerformanceCounter::StartCounter()
{
  struct timeval tv;

  gettimeofday(&tv,NULL);

  startCountSec = tv.tv_sec;
  startCountMicroSec = tv.tv_usec;
}

inline void PerformanceCounter::StopCounter()
{
  struct timeval tv;

  gettimeofday(&tv,NULL);

  stopCountSec = tv.tv_sec;
  stopCountMicroSec = tv.tv_usec;
}


inline double PerformanceCounter::GetElapsedTime()
{
  float elapsedTime = 1.0 * (stopCountSec-startCountSec) + 1E-6 * (stopCountMicroSec - startCountMicroSec);
  return elapsedTime;
}

#endif

#ifdef WIN32

/**************** WINDOWS COUNTER *******************/

#include <windows.h>

class PerformanceCounter
{
public:

  PerformanceCounter()
  {
    // reset the counter frequency
    QueryPerformanceFrequency(&timerFrequency);
    StartCounter();
  }

  void StartCounter(); // call this before your code block
  void StopCounter(); // call this after your code block

  // read elapsed time (units are seconds, accuracy is up to microseconds)
  double GetElapsedTime();

protected:
  LARGE_INTEGER timerFrequency;
  LARGE_INTEGER startCount,stopCount;
};

inline void PerformanceCounter::StartCounter()
{
  QueryPerformanceCounter(&startCount);
}

inline void PerformanceCounter::StopCounter()
{
  QueryPerformanceCounter(&stopCount);
}


inline double PerformanceCounter::GetElapsedTime()
{
  return ((double)(stopCount.QuadPart - startCount.QuadPart))
    / ((double)timerFrequency.QuadPart);
}

#endif
#endif
/*
skeleton.h

Definition of the skeleton.

Written by Jehee Lee

Revision 1 - Steve Lin, Jan. 14, 2002
Revision 2 - Alla and Kiran, Jan 18, 2002
Revision 3 - Jernej Barbic and Yili Zhao, Feb, 2012
*/

#ifndef _VECTOR_H
#define _VECTOR_H


class vector
{
  // negation
  friend vector    operator-( vector const& );

  // addtion
  friend vector    operator+( vector const&, vector const& );

  // subtraction
  friend vector    operator-( vector const&, vector const& );

  // dot product
  friend double    operator%( vector const&, vector const& );

  // cross product
  friend vector    operator*( vector const&, vector const& );

  // scalar Multiplication
  friend vector    operator*( vector const&, double );

  // scalar Division
  friend vector    operator/( vector const&, double );


  friend double    len( vector const& );
  friend vector	normalize( vector const& );

  friend double       angle( vector const&, vector const& );

  // member functions
public:
  // constructors
  vector() {}
  vector( double x, double y, double z ) { p[0]=x; p[1]=y; p[2]=z; }
  vector( double a[3] ) { p[0]=a[0]; p[1]=a[1]; p[2]=a[2]; }
  ~vector() {};

  // inquiry functions
  double& operator[](int i) { return p[i];}

  double x() const { return p[0]; };
  double y() const { return p[1]; };
  double z() const { return p[2]; };
  void   getValue( double d[3] ) { d[0]=p[0]; d[1]=p[1]; d[2]=p[2]; }
  void   setValue( double d[3] ) { p[0]=d[0]; p[1]=d[1]; p[2]=d[2]; }

  double getValue( int n ) const { return p[n]; }
  vector setValue( double x, double y, double z )
  { p[0]=x, p[1]=y, p[2]=z; return *this; }
  double setValue( int n, double x )
  { return p[n]=x; }

  double length() const;

  // change functions
  void set_x( double x ) { p[0]=x; };
  void set_y( double x ) { p[1]=x; };
  void set_z( double x ) { p[2]=x; };

  //data members
  double p[3]; //X, Y, Z components of the vector
};

#endif
