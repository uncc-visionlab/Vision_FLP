#!/bin/bash
export LD_LIBRARY_PATH=/usr/lib
source /opt/ros/kinetic/setup.bash
source ../../../devel/setup.bash
export ROS_NAMESPACE=/camera
rosrun ptcloud_engine py_rgbd_listener.py _filename:=$1