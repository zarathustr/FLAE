# FLAE
The Fast Linear Attitude Estimator (FLAE) proposed by Jin Wu et al. in Feb-2016 during the undergraduate study of Jin Wu. This estimator is aimed to compute the attitude of rigid body (like satellite, vehicles, smart devices and etc.) using 3D vector measurements including star vectors, sun vectors, lunar vectors, geomagnetic vectors, gravitational vectors, visual-tracking points and etc.

# Citation
Wu, J., Zhou, Z., Gao, B., Li, R., Cheng, Y., & Fourati, H. (2018). Fast Linear Quaternion Attitude Estimator Using Vector Observations. IEEE Transactions on Automation Science and Engineering, 15(1), 307–319. https://doi.org/10.1109/TASE.2017.2699221 (**ESI Highly Cited Paper by ISI Web of Science in 2017-2018 year round**).

# Usage
In MATLAB software, type:

**$ q = FLAE(body_vectors, ref_vectors, weights, 'symbolic');**

in which **body_vectors, ref_vectors** are stacked body and reference vector matrices with 3xn dimensions; **weights** is an 1xn or nx1 vector with all positive numbers and their sum equals to 1; The option **'symbolic'** points to the solving method using symbolic solutions to the quatic polynomial while it can be replaced by **'eig'** and **'newton'** as well which correspond to the eigen-decomposition of MATLAB and Newton iteration respectively. The final result is in the form of an unitary quaternion **q** that can be later cast into rotation matrix or Euler angles.

# Acknowledgement
The authors would like to thank Dr. Yaguang Yang from NRC and Dr. F. Landis Markley from NASA Goddard Space Flight Center for their comments and assistance in improving the architecture of this algorithm.
