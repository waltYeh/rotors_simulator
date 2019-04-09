#!/usr/bin/env python
# rosrun rotors_gazebo waypoint_pub_file.py $(find rotors_gazebo)/resource/example_waypoints.txt
from trajectory_msgs.msg import MultiDOFJointTrajectory, MultiDOFJointTrajectoryPoint
from geometry_msgs.msg import Vector3
from geometry_msgs.msg import Twist
from geometry_msgs.msg import Transform, Quaternion
import mav_msgs.msg as mav_msgs
from sensor_msgs.msg import Imu
import tf
import rospy
import sys
from rotors_gazebo.srv import *
DEG_2_RAD = 0.0174532925
kNanoSecondsInSecond = 1000000000
sim_running = False
def callback(data):
	global sim_running
	if sim_running == False:
		rospy.loginfo("sim_running")
	sim_running = True
def traj_gen_client():
	rospy.wait_for_service('traj_gen')
	try:
		traj_gen = rospy.ServiceProxy('traj_gen', traj_gen_service)
		resp1 = traj_gen(True)
		return resp1.trj
	except rospy.ServiceException, e:
		print "Service call failed: %s"%e
if __name__=='__main__':
	rospy.init_node('waypoint_publisher', anonymous=True)

	rospy.Subscriber("/hummingbird/imu", Imu, callback)
	pub = rospy.Publisher('/hummingbird/command/trajectory', MultiDOFJointTrajectory,queue_size=10)
	rospy.loginfo("Wait for simulation to become ready...")
	while (not sim_running ) and (not rospy.is_shutdown()):
		rospy.sleep(1.)

	rospy.loginfo("ok")
	rospy.sleep(5.)
	rospy.loginfo("Start computing waypoints.")

#	msg=MultiDOFJointTrajectory()
	# msg.header.frame_id =''
	# msg.header.stamp = rospy.Time.now()
	# msg.joint_names = ["base_link"]
	msg = traj_gen_client()
	rospy.loginfo("Start publishing waypoints.")
	pub.publish(msg)