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
DEG_2_RAD = 0.0174532925
kNanoSecondsInSecond = 1000000000
sim_running = False
def callback(data):
	global sim_running
	if sim_running == False:
		rospy.loginfo("sim_running")
	sim_running = True
class WaypointWithTime:
	def __init__(self):
		self.yaw = 0.
		self.waiting_time = 0.
		self.position = Vector3()
		self.position.x = 0.0
		self.position.y = 0.0
		self.position.z = 0.0
	def __init__(self,t,x,y,z,_yaw):
		self.yaw = _yaw
		self.waiting_time = t
		self.position = Vector3()
		self.position.x = x
		self.position.y = y
		self.position.z = z
if __name__=='__main__':
	# print(sys.argv[0])
	# print(len(sys.argv))
	# print(str(sys.argv))
	rospy.init_node('waypoint_publisher', anonymous=True)
	waypoints = []
	file = open(sys.argv[1],"r")
	for linestr in file.readlines():
		linelist=[]
		linelist+=linestr.split()
		time = float(linelist[0])
		x = float(linelist[1])
		y = float(linelist[2])
		z = float(linelist[3])
		yaw = float(linelist[4])
		waypoints += [WaypointWithTime(time,x,y,z,yaw * DEG_2_RAD)]
	file.close()
	rospy.loginfo("Read %d waypoints.", len(waypoints))
	rospy.Subscriber("/hummingbird/imu", Imu, callback)
	pub = rospy.Publisher('/hummingbird/command/trajectory', MultiDOFJointTrajectory,queue_size=10)
	rospy.loginfo("Wait for simulation to become ready...")
	while (not sim_running ) and (not rospy.is_shutdown()):
		rospy.sleep(1.)

	rospy.loginfo("ok")
	rospy.sleep(5.)
	rospy.loginfo("Start publishing waypoints.")
	msg=MultiDOFJointTrajectory()
	msg.header.frame_id =''
	msg.header.stamp = rospy.Time.now()
	msg.joint_names = ["base_link"]
	time_from_start = rospy.Duration.from_sec(0)
	for i in range(len(waypoints)):
		transforms = Transform()
		transforms.translation.x = waypoints[i].position.x
		transforms.translation.y = waypoints[i].position.y
		transforms.translation.z = waypoints[i].position.z

		quaternion = tf.transformations.quaternion_from_euler(0, 0, waypoints[i].yaw)

		transforms.rotation=Quaternion(quaternion[0],quaternion[1],quaternion[2],quaternion[3])

		velocities =Twist()
		accelerations=Twist()
		time_from_start += rospy.Duration.from_sec(waypoints[i].waiting_time)
		point = MultiDOFJointTrajectoryPoint([transforms],[velocities],[accelerations],time_from_start)
		
		msg.points.append(point)
	pub.publish(msg)