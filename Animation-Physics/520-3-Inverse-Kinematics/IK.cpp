#include "IK.h"
#include "FK.h"
#include "minivectorTemplate.h"
#include <Eigen/Dense>
#include <adolc/adolc.h>
#include <cassert>
#if defined(_WIN32) || defined(WIN32)
  #ifndef _USE_MATH_DEFINES
    #define _USE_MATH_DEFINES
  #endif
#endif
#include <math.h>
using namespace std;

// CSCI 520 Computer Animation and Simulation
// Jernej Barbic and Yijing Li

namespace
{

// Converts degrees to radians.
template<typename real>
inline real deg2rad(real deg) { return deg * M_PI / 180.0; }

template<typename real>
Mat3<real> Euler2Rotation(const real angle[3], RotateOrder order)
{
  Mat3<real> RX = Mat3<real>::getElementRotationMatrix(0, deg2rad(angle[0]));
  Mat3<real> RY = Mat3<real>::getElementRotationMatrix(1, deg2rad(angle[1]));
  Mat3<real> RZ = Mat3<real>::getElementRotationMatrix(2, deg2rad(angle[2]));

  switch(order)
  {
    case RotateOrder::XYZ:
      return RZ * RY * RX;
    case RotateOrder::YZX:
      return RX * RZ * RY;
    case RotateOrder::ZXY:
      return RY * RX * RZ;
    case RotateOrder::XZY:
      return RY * RZ * RX;
    case RotateOrder::YXZ:
      return RZ * RX * RY;
    case RotateOrder::ZYX:
      return RX * RY * RZ;
  }
  assert(0);
}

// Performs forward kinematics, using the provided "fk" class.
// This is the function whose Jacobian matrix will be computed using adolc.
// numIKJoints and IKJointIDs specify which joints serve as handles for IK:
//   IKJointIDs is an array of integers of length "numIKJoints"
// Input: numIKJoints, IKJointIDs, fk, eulerAngles (of all joints)
// Output: handlePositions (world-coordinate positions of all the IK joints; length is 3 * numIKJoints)
template<typename real>
void forwardKinematicsFunction(
    int numIKJoints, const int* IKJointIDs, const FK& fk,
    const std::vector<real>& eulerAngles, std::vector<real>& handlePositions)
{
    // Students should implement this.
    // The implementation of this function is very similar to function computeLocalAndGlobalTransforms in the FK class.
    // The recommended approach is to first implement FK::computeLocalAndGlobalTransforms.
    // Then, implement the same algorithm into this function. To do so,
    // you can use fk.getJointUpdateOrder(), fk.getJointRestTranslation(), and fk.getJointRotateOrder() functions.
    // Also useful is the multiplyAffineTransform4ds function in minivectorTemplate.h .
    // It would be in principle possible to unify this "forwardKinematicsFunction" and FK::computeLocalAndGlobalTransforms(),
    // so that code is only written once. We considered this; but it is actually not easily doable.
    // If you find a good approach, feel free to document it in the README file, for extra credit.
    int numJoints = fk.getNumJoints();

    std::vector<Mat3<real>> localR(numJoints), globalR(numJoints);
    std::vector<Vec3<real>> localT(numJoints), globalT(numJoints);

    // local transforms
    for (int i = 0; i < numJoints; ++i)
    {
        // joint rotation
        real angle[3] = {
            eulerAngles[3 * i + 0],
            eulerAngles[3 * i + 1],
            eulerAngles[3 * i + 2]
        };
        Mat3<real> R = Euler2Rotation(angle, fk.getJointRotateOrder(i));

        // joint orientation rotation
        Vec3d jOrient = fk.getJointOrient(i);
        real orientAngle[3] = {
            static_cast<real>(jOrient[0]),
            static_cast<real>(jOrient[1]),
            static_cast<real>(jOrient[2])
        };
        Mat3<real> R_orient = Euler2Rotation(orientAngle, fk.getJointRotateOrder(i)); // assuming same order

        // local rotation
        localR[i] = R_orient * R;

        // local translation
        Vec3d restT = fk.getJointRestTranslation(i);
        localT[i] = Vec3<real>(restT[0], restT[1], restT[2]);
    }

    // global transforms
    for (int i = 0; i < numJoints; ++i)
    {
        int jointID = fk.getJointUpdateOrder(i);
        int parentID = fk.getJointParent(jointID);

        if (parentID < 0) {
            globalR[jointID] = localR[jointID];
            globalT[jointID] = localT[jointID];
        }
        else {
            multiplyAffineTransform4ds(
                globalR[parentID], globalT[parentID],
                localR[jointID], localT[jointID],
                globalR[jointID], globalT[jointID]);
        }
    }

    //  handle positions
    handlePositions.resize(numIKJoints * 3);
    for (int i = 0; i < numIKJoints; ++i)
    {
        int jointID = IKJointIDs[i];
        handlePositions[3 * i + 0] = globalT[jointID][0];
        handlePositions[3 * i + 1] = globalT[jointID][1];
        handlePositions[3 * i + 2] = globalT[jointID][2];
    }
}
} // end anonymous namespaces

IK::IK(int numIKJoints, const int * IKJointIDs, FK * inputFK, int adolc_tagID)
{
  this->numIKJoints = numIKJoints;
  this->IKJointIDs = IKJointIDs;
  this->fk = inputFK;
  this->adolc_tagID = adolc_tagID;

  FKInputDim = fk->getNumJoints() * 3;
  FKOutputDim = numIKJoints * 3;

  train_adolc();
}

void IK::train_adolc()
{
  // Students should implement this.
  // Here, you should setup adol_c:
  //   Define adol_c inputs and outputs. 
  //   Use the "forwardKinematicsFunction" as the function that will be computed by adol_c.
  //   This will later make it possible for you to compute the gradient of this function in IK::doIK
  //   (in other words, compute the "Jacobian matrix" J).
  // See ADOLCExample.cpp .
    
    trace_on(adolc_tagID); // start tracking

    vector<adouble> eulerAngles(FKInputDim);
    for (int i = 0; i < FKInputDim; ++i)
        eulerAngles[i] <<= 0.0; // Input to ADOL-C

    vector<adouble> handlePositions(FKOutputDim);
    forwardKinematicsFunction(numIKJoints, IKJointIDs, *fk, eulerAngles, handlePositions); // Call templated FK

    vector<double> output(FKOutputDim);
    for (int i = 0; i < FKOutputDim; ++i)
        handlePositions[i] >>= output[i]; // Output 

    trace_off(); // stop

}

void IK::doIK(const Vec3d * targetHandlePositions, Vec3d * jointEulerAngles)
{
  // You may find the following helpful:
  int numJoints = fk->getNumJoints(); // Note that is NOT the same as numIKJoints!

  // Students should implement this.
  // Use adolc to evalute the forwardKinematicsFunction and its gradient (Jacobian). It was trained in train_adolc().
  // Specifically, use ::function, and ::jacobian .
  // See ADOLCExample.cpp .
  //
  // Use it implement the Tikhonov IK method (or the pseudoinverse method for extra credit).
  // Note that at entry, "jointEulerAngles" contains the input Euler angles. 
  // Upon exit, jointEulerAngles should contain the new Euler angles.


  vector<double> x(FKInputDim);
  for (int i = 0; i < numJoints; ++i)
  {
      x[i * 3 + 0] = jointEulerAngles[i][0];
      x[i * 3 + 1] = jointEulerAngles[i][1];
      x[i * 3 + 2] = jointEulerAngles[i][2];
  }


  vector<double> y(FKOutputDim);
  ::function(adolc_tagID, FKOutputDim, FKInputDim, x.data(), y.data());

  // create target
  vector<double> target(FKOutputDim);
  for (int i = 0; i < numIKJoints; ++i)
  {
      target[i * 3 + 0] = targetHandlePositions[i][0];
      target[i * 3 + 1] = targetHandlePositions[i][1];
      target[i * 3 + 2] = targetHandlePositions[i][2];
  }

  // diff = target - y
  vector<double> diff(FKOutputDim);
  for (int i = 0; i < FKOutputDim; ++i)
      diff[i] = target[i] - y[i];

  // jacobian
  vector<double> jacobian(FKInputDim * FKOutputDim);
  vector<double*> jacobianRows(FKOutputDim);
  for (int i = 0; i < FKOutputDim; ++i)
      jacobianRows[i] = &jacobian[i * FKInputDim];
  ::jacobian(adolc_tagID, FKOutputDim, FKInputDim, x.data(), jacobianRows.data());

  Eigen::MatrixXd J(FKOutputDim, FKInputDim);
  for (int i = 0; i < FKOutputDim; ++i)
      for (int j = 0; j < FKInputDim; ++j)
          J(i, j) = jacobian[i * FKInputDim + j];

  Eigen::VectorXd b(FKOutputDim);
  for (int i = 0; i < FKOutputDim; ++i)
      b[i] = diff[i];

  double lambda = 1e-3;
  Eigen::MatrixXd JTJ = J.transpose() * J;
  Eigen::VectorXd JTb = J.transpose() * b;
  Eigen::VectorXd dx = (JTJ + lambda * Eigen::MatrixXd::Identity(FKInputDim, FKInputDim)).ldlt().solve(JTb);

  // jointEulerAngles
  for (int i = 0; i < numJoints; ++i)
  {
      jointEulerAngles[i][0] += dx[i * 3 + 0];
      jointEulerAngles[i][1] += dx[i * 3 + 1];
      jointEulerAngles[i][2] += dx[i * 3 + 2];
  }

}

