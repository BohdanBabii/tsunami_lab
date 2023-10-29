.. _ch:Task_1.1:

Report Week 1
==============
To implement the f-wave solver we used a class called F_Wave the only public method of F_Wave is called netUpdates, which calculates all netUpdates as required in the task. To calculate the netUpdates we need to compute each wave first, that again requires the wave speeds and wave strengths to be available. For which we wrote the private functions waveSpeeds and waveStrengths. Lastly we added the heightAvg and the particleVelocityAvg function to be able to calculate the wave speeds. 

Both heightAvg, which calculates average height of the left and the right cell, and particleVelocityAvg, computing the average particle velocity of both cells, are implemented as extra functions to increase readability.
This chapter implements and tests the most basic part of our project: The *f-wave* solver for the one-dimensional shallow water equations.
The shallow water equations are a system of nonlinear hyperbolic conservations laws with an optional source term:

.. math:: h(q_l, q_r) = \frac{1}{2}(h_l+h_r), u(q_l, q_r) = \frac{u_l\sqrt{h_l}+u_r\sqrt{h_r}}{\sqrt{h_l}+\sqrt{h_r}}.

The waveSpeeds function computes the roe eigenvalues (wave speeds of both waves) using the average height, average particle velocity and the gravity constant g.

.. math:: \lambda_{1, 2}=u(q_l, q_r)\mp\sqrt{g\cdot h(q_l, q_r)}

The waveStrengths function computes the wave strengths of both waves.

.. math:: \begin{bmatrix}\alpha_1 \\ \alpha_2 \end{bmatrix} = \begin{bmatrix}1 & 1\\ \lambda_1 & \lambda_2\end{bmatrix}^{-1}\Delta f 

To increase the readability we used the Variable inv_det that stores the inverse determinant.


.. math:: \frac{1}{\lambda_2-\lambda_1}\begin{bmatrix}\lambda_2 & -1\\ -\lambda_1 & 1\end{bmatrix} = \begin{bmatrix}\lambda_2\cdot inv\_det & -inv\_det\\ -\lambda_1\cdot inv\_det & inv\_det\end{bmatrix}\qquad

When we wrote out the matrix-vector multiplication we got formulas for each wave strength, which can be calculated easily by a computer.

.. math:: \alpha_1 = \lambda_2\cdot inv\_det\cdot (h_l- h_r) - inv\_det\cdot(hu_l-hu_r),\\ \alpha_2 = inv\_det\cdot(hu_l-hu_r)-\lambda_1\cdot inv\_det\cdot(h_l-h_r)

If the wave strength is greater than 0 it belongs to a wave is right-going which influences the right cell. And the other way around. We expect one wave to be right-going and on to be left-going. This way the left net update equals the left-going wave and the right net update is set to be the right going wave.

.. math:: netUpdate_{Left}= A^{-}\Delta Q = \begin{cases}Z_1\qquad\lambda_1<0\\ Z_2\qquad\lambda_2<0\end{cases} \\ netUpdate_{Right}= A^{+}\Delta Q = \begin{cases}Z_1\qquad\lambda_1>0\\ Z_2\qquad\lambda_2>0\end{cases} 


.. _ch:code:

Code
---------------
.. code-block:: cpp


	/**
	 * @author Bohdan Babii, Phillip Rothenbeck
	 *
	 * @section DESCRIPTION
	 * F Wave solver for the one-dimensional shallow water equations.
	**/

	#ifndef TSUNAMI_LAB_SOLVERS_F_WAFE
	#define TSUNAMI_LAB_SOLVERS_F_WAFE

	#include "../constants.h"

	namespace tsunami_lab {
		namespace solvers {
			class F_Wave;
		}
	}

	class tsunami_lab::solvers::F_Wave {
		private:

			//! square root of gravity
			static t_real constexpr c_sqrt_g = 3.131557121;

			/**
			 * Computes the average wave height.
			 *
			 * @param i_hL height of the left side.
			 * @param i_hR height of the right side.
			 * @param o_hight will be set to the average speed.
			**/

			static void heightAvg( t_real 	i_hL,
								   t_real 	i_hR,
								   t_real & o_height);

	        /**
			 * Computes the average particle_Velocity
			 *
	         * @param i_hL height of the left side.
	         * @param i_hR height of the right side.
	         * @param i_huL momentum of the left side.
	         * @param i_huR momentum of the right side.
	         * @param o_velocity will be set to the average velocity.
	        **/

			static void particleVelocityAvg( t_real  i_hL,
											 t_real  i_hR,
											 t_real  i_uL,
											 t_real	 i_uR,
											 t_real & o_velocity);
	        /**
			 * Computes the wave speeds.
			 *
	         * @param i_hL height of the left side.
	         * @param i_hR height of the right side.
	         * @param i_huL momentum of the left side.
	         * @param i_huR momentum of the right side.
	         * @param o_speed_left will be set to the speed of the wave propagating to the left.
	         * @param o_speed_right will be set to the speed of the wave propagating to the right.
	        **/

			static void waveSpeeds(	t_real   i_hL,
									t_real   i_hR,
									t_real   i_uL,
									t_real   i_uR,
									t_real & o_wafeSpeedL,
									t_real & o_wafeSpeedR);

	        /**
			 * Computes the wave strengths
			 * 
	         * @param i_hL height of the left side.
	         * @param i_hR height of the right side.
	         * @param i_huL momentum of the left side.
	         * @param i_huR momentum of the right side.
	         * @param o_waveSpeeds will be set to the strength of the wave propagation to the left.
			 * @param o_wafeSpeeds will be set to the strength of the wave propagation to the right.
	        **/

			static void waveStrengths( t_real   i_hL,
									   t_real   i_hR,
									   t_real   i_huL,
									   t_real   i_huR,
									   t_real   i_waveSpeedL,
	                               	   t_real   i_waveSpeedR,
									   t_real & o_strengthL,
									   t_real & o_strengthR);

		public:
	        /**
			 * Computes the net-updates.
			 *
	         * @param i_hL height of the left side.
	         * @param i_hR height of the right side.
	         * @param i_huL momentum of the left side.
	         * @param i_huR momentum of the right side.
	         * @param o_netUpdateL will be set to the net-updates for the left side; 0: hight 1: momentum.
			 * @param o_netUpdateR will be set to the net-updates for the right side; 0: hight, 1: momentum. 
	        **/

			static void netUpdates( t_real i_hL,
	                            	t_real i_hR,
	                            	t_real i_huL,
	                           		t_real i_huR,
	                            	t_real o_netUpdateL[2],
	                            	t_real o_netUpdateR[2] );
	};							
	#endif
