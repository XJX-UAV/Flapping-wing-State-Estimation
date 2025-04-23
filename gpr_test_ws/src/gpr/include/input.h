#include <ros/ros.h>
#include <Eigen/Dense>

#include <sensor_msgs/Imu.h>
#include <geometry_msgs/TwistStamped.h>
#include <sensor_msgs/MagneticField.h>
#include <nav_msgs/Odometry.h>
#include <mavros_msgs/RCIn.h>

class Imu_Data_t
{
public:
  Eigen::Quaterniond q;
  Eigen::Vector3d w;
  Eigen::Vector3d a;

  sensor_msgs::Imu msg;
  ros::Time rcv_stamp;

  // Imu_Data_t();
  void feed(sensor_msgs::ImuConstPtr pMsg);
};

class Odom_Data_t
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  Eigen::Vector3d p;
  Eigen::Vector3d v;
  Eigen::Quaterniond q;

  nav_msgs::Odometry msg;
  ros::Time rcv_stamp;
  bool recv_new_msg;

  Odom_Data_t();
  void feed(nav_msgs::OdometryConstPtr pMsg);
};

class RC_Data_t
{
public:
  bool est_trigger;

  mavros_msgs::RCIn msg;
  ros::Time rcv_stamp;

  void feed(mavros_msgs::RCInConstPtr pMsg);
};

class GPS_VEL_Data_t
{
public:
  Eigen::Vector3d gps_vel;

  geometry_msgs::TwistStamped msg;
  ros::Time rcv_stamp;

  void feed(geometry_msgs::TwistStampedConstPtr pMsg);
};

class Mag_Data_t
{
public:
  Eigen::Vector3d mag;

  sensor_msgs::MagneticField msg;
  ros::Time rcv_stamp;

  void feed(sensor_msgs::MagneticFieldConstPtr pMsg);
};