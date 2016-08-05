// FCF_test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

std::ifstream input1("flat_motion2.txt");
std::ifstream input2("flat_motion2.txt");
std::ofstream output("output.txt");


float invSqrt(float x)
{
	float xhalf = 0.5f * x;
	int i = *(int*)&x; // get bits for floating value
	i = 0x5f375a86 - (i >> 1); // gives initial guess
	x = *(float*)&i; // convert bits back to float
	x = x * (1.5f - xhalf*x*x); // Newton step
	return x;
}


void acc_mag(const float M[6], const float q_a[4], float q[4])
{
	int k;
	static const signed char iv0[4] = { 1, 0, 0, 0 };

	float ab_idx_0;
	float ab_idx_1;
	float ab_idx_2;
	float ab_idx_3;
	float h_idx_1;
	float h_idx_2;
	float h_idx_3;
	float mn;
	float b1xr1[3];
	float b2xr2[3];
	float b_M[3];
	float bxrx;
	float c_M[3];
	float bxxrx[3];
	float bx_rx[3];
	float b_bxxrx;
	float alpha;
	float beta;
	float b_gamma;
	float norma;
	for (k = 0; k < 4; k++) {
		q[k] = iv0[k];
	}

	/* QUATERN2ROTMAT Converts a quaternion to its conjugate */
	/*  */
	/*    qConj = quaternConj(q) */
	/*  */
	/*    Converts a quaternion to its conjugate. */
	/*  */
	/*    For more information see: */
	/*    http://www.x-io.co.uk/node/8#quaternions */
	/*  */
	/* 	Date          Author          Notes */
	/* 	27/09/2011    SOH Madgwick    Initial release */
	/* QUATERNPROD Calculates the quaternion product */
	/*  */
	/*    ab = quaternProd(a, b) */
	/*  */
	/*    Calculates the quaternion product of quaternion a and b. */
	/*  */
	/*    For more information see: */
	/*    http://www.x-io.co.uk/node/8#quaternions */
	/*  */
	/* 	Date          Author          Notes */
	/* 	27/09/2011    SOH Madgwick    Initial release */
	ab_idx_0 = ((0.0f * q_a[0] - M[3] * -q_a[1]) - M[4] * -q_a[2]) - M[5] * -q_a[3];
	ab_idx_1 = ((0.0f * -q_a[1] + M[3] * q_a[0]) + M[4] * -q_a[3]) - M[5] * -q_a[2];
	ab_idx_2 = ((0.0f * -q_a[2] - M[3] * -q_a[3]) + M[4] * q_a[0]) + M[5] * -q_a[1];
	ab_idx_3 = ((0.0f * -q_a[3] + M[3] * -q_a[2]) - M[4] * -q_a[1]) + M[5] * q_a[0];

	/* QUATERNPROD Calculates the quaternion product */
	/*  */
	/*    ab = quaternProd(a, b) */
	/*  */
	/*    Calculates the quaternion product of quaternion a and b. */
	/*  */
	/*    For more information see: */
	/*    http://www.x-io.co.uk/node/8#quaternions */
	/*  */
	/* 	Date          Author          Notes */
	/* 	27/09/2011    SOH Madgwick    Initial release */
	h_idx_1 = ((q_a[0] * ab_idx_1 + q_a[1] * ab_idx_0) + q_a[2] * ab_idx_3) - q_a
		[3] * ab_idx_2;
	h_idx_2 = ((q_a[0] * ab_idx_2 - q_a[1] * ab_idx_3) + q_a[2] * ab_idx_0) + q_a
		[3] * ab_idx_1;
	h_idx_3 = ((q_a[0] * ab_idx_3 + q_a[1] * ab_idx_2) - q_a[2] * ab_idx_1) + q_a
		[3] * ab_idx_0;
	mn = sqrt(h_idx_1 * h_idx_1 + h_idx_2 * h_idx_2);
	b1xr1[0] = M[1];
	b1xr1[1] = -M[0];
	b1xr1[2] = 0.0f;
	b2xr2[0] = M[4] * h_idx_3;
	b2xr2[1] = -M[3] * h_idx_3 + M[5] * mn;
	b2xr2[2] = -M[4] * mn;
	b_M[0] = M[1] * M[3] - M[0] * M[4];
	b_M[1] = 0.0f;
	b_M[2] = -M[2] * M[4] + M[1] * M[5];
	bxrx = (M[2] * M[3] - M[0] * M[5]) / mn;
	c_M[0] = -M[2] * M[4] + M[1] * M[5];
	c_M[1] = (M[2] * M[3] - M[0] * M[5]) + mn;
	c_M[2] = -M[1] * M[3] + M[0] * M[4];
	ab_idx_0 = 0.0f;
	for (k = 0; k < 3; k++) {
		b_bxxrx = b_M[k] / mn;
		bx_rx[k] = c_M[k] / mn;
		ab_idx_0 += b_bxxrx * (0.5f * b1xr1[k] + 0.5f * b2xr2[k]);
		bxxrx[k] = b_bxxrx;
		b1xr1[k] = 0.5f * b1xr1[k] + 0.5f * b2xr2[k];
	}

	alpha = (1.0f + bxrx) * (0.5f * M[2] + 0.5f * (M[5] * h_idx_3 + M[3] * mn)) +
		ab_idx_0;
	beta = 0.0f;
	for (k = 0; k < 3; k++) {
		beta += bx_rx[k] * b1xr1[k];
	}

	b_gamma = sqrt(alpha * alpha + beta * beta);
	if (alpha >= 0.0) {
		norma = 2.0f / invSqrt(b_gamma * (b_gamma + alpha) * (1.0f + bxrx));
		ab_idx_0 = b_gamma + alpha;
		q[0] = (b_gamma + alpha) * (1.0f + bxrx) / norma;
		for (k = 0; k < 3; k++) {
			q[k + 1] = (ab_idx_0 * bxxrx[k] + beta * bx_rx[k]) / norma;
		}
	}

	if (alpha < 0.0) {
		norma = 2.0f / invSqrt(b_gamma * (b_gamma - alpha) * (1.0f + bxrx));
		ab_idx_0 = b_gamma - alpha;
		q[0] = beta * (1.0f + bxrx) / norma;
		for (k = 0; k < 3; k++) {
			q[k + 1] = (beta * bxxrx[k] + ab_idx_0 * bx_rx[k]) / norma;
		}
	}
}





void FCF_gyro_acc_mag(float *Gyroscope, float *Accelerometer, float *Magnetometer, float *q_,
	float gain_a, float gain_m, float std_norm_m, float threshold_m, float dt,float *q_est)
{

	float norm_a, norm_m,q_acc[4],q_gravity_mag[4],norm_q;

	norm_a = invSqrt(Accelerometer[0] * Accelerometer[0] + Accelerometer[1] * Accelerometer[1] + Accelerometer[2] * Accelerometer[2]);
	Accelerometer[0] *= norm_a;
	Accelerometer[1] *= norm_a;
	Accelerometer[2] *= norm_a;

	norm_m = invSqrt(Magnetometer[0] * Magnetometer[0] + Magnetometer[1] * Magnetometer[1] + Magnetometer[2] * Magnetometer[2]);
	Magnetometer[0] *= norm_m;
	Magnetometer[1] *= norm_m;
	Magnetometer[2] *= norm_m;

	float P[4][4] = {
		{ Accelerometer[2] + 1.0f	, Accelerometer[1]			, -Accelerometer[0]			, 0							},
		{ Accelerometer[1]			, -Accelerometer[2] + 1.0f	, 0							, Accelerometer[0]			},
		{ -Accelerometer[0]			, 0							, -Accelerometer[2] + 1.0f	, Accelerometer[1]			},
		{ 0							, Accelerometer[0]			, Accelerometer[1]			, Accelerometer[2] + 1.0f	}
	};

	float omega[4][4] = {
		{ 0				, -Gyroscope[0]	, -Gyroscope[1]	, -Gyroscope[2]    },
		{ Gyroscope[0]	, 0				, Gyroscope[2]	, -Gyroscope[1]    },
		{ Gyroscope[1]	, -Gyroscope[2]	, 0				, Gyroscope[0]	   },
		{ Gyroscope[2]	, Gyroscope[1]	, -Gyroscope[0]	,				 0 }
	};

	if (fabs(1.0f/norm_a - 1.0f) < 0.2f)
	{

		float q_0[4], q_1[4];

		q_0[0] = (0.5f*dt*omega[0][0] + 1.0f) * q_[0] + (0.5f*dt*omega[0][1]) * q_[1] + (0.5f*dt*omega[0][2]) * q_[2] + (0.5f*dt*omega[0][3]) * q_[3];
		q_0[1] = (0.5f*dt*omega[1][0]) * q_[0] + (0.5f*dt*omega[1][1] + 1.0f) * q_[1] + (0.5f*dt*omega[1][2]) * q_[2] + (0.5f*dt*omega[1][3]) * q_[3];
		q_0[2] = (0.5f*dt*omega[2][0]) * q_[0] + (0.5f*dt*omega[2][1]) * q_[1] + (0.5f*dt*omega[2][2] + 1.0f) * q_[2] + (0.5f*dt*omega[2][3]) * q_[3];
		q_0[3] = (0.5f*dt*omega[3][0]) * q_[0] + (0.5f*dt*omega[3][1]) * q_[1] + (0.5f*dt*omega[3][2]) * q_[2] + (0.5f*dt*omega[3][3] + 1.0f) * q_[3];

		q_0[0] *= (1.0f - gain_a);
		q_0[1] *= (1.0f - gain_a);
		q_0[2] *= (1.0f - gain_a);
		q_0[3] *= (1.0f - gain_a);

		q_1[0] = gain_a*0.5f*P[0][0] * q_[0] + gain_a*0.5f*P[0][1] * q_[1] + gain_a*0.5f*P[0][2] * q_[2] + gain_a*0.5f*P[0][3] * q_[3];
		q_1[1] = gain_a*0.5f*P[1][0] * q_[0] + gain_a*0.5f*P[1][1] * q_[1] + gain_a*0.5f*P[1][2] * q_[2] + gain_a*0.5f*P[1][3] * q_[3];
		q_1[2] = gain_a*0.5f*P[2][0] * q_[0] + gain_a*0.5f*P[2][1] * q_[1] + gain_a*0.5f*P[2][2] * q_[2] + gain_a*0.5f*P[2][3] * q_[3];
		q_1[3] = gain_a*0.5f*P[3][0] * q_[0] + gain_a*0.5f*P[3][1] * q_[1] + gain_a*0.5f*P[3][2] * q_[2] + gain_a*0.5f*P[3][3] * q_[3];


		q_acc[0] = q_0[0] + q_1[0];
		q_acc[1] = q_0[1] + q_1[1];
		q_acc[2] = q_0[2] + q_1[2];
		q_acc[3] = q_0[3] + q_1[3];

	}
	else
	{
		q_acc[0] = (0.5f*dt*omega[0][0] + 1.0f) * q_[0] + (0.5f*dt*omega[0][1]) * q_[1] + (0.5f*dt*omega[0][2]) * q_[2] + (0.5f*dt*omega[0][3]) * q_[3];
		q_acc[1] = (0.5f*dt*omega[1][0]) * q_[0] + (0.5f*dt*omega[1][1] + 1.0f) * q_[1] + (0.5f*dt*omega[1][2]) * q_[2] + (0.5f*dt*omega[1][3]) * q_[3];
		q_acc[2] = (0.5f*dt*omega[2][0]) * q_[0] + (0.5f*dt*omega[2][1]) * q_[1] + (0.5f*dt*omega[2][2] + 1.0f) * q_[2] + (0.5f*dt*omega[2][3]) * q_[3];
		q_acc[3] = (0.5f*dt*omega[3][0]) * q_[0] + (0.5f*dt*omega[3][1]) * q_[1] + (0.5f*dt*omega[3][2]) * q_[2] + (0.5f*dt*omega[3][3] + 1.0f) * q_[3];
	}

	norm_q = invSqrt(q_acc[0] * q_acc[0] + q_acc[1] * q_acc[1] + q_acc[2] * q_acc[2] + q_acc[3] * q_acc[3]);
	q_acc[0] *= norm_q;
	q_acc[1] *= norm_q;
	q_acc[2] *= norm_q;
	q_acc[3] *= norm_q;

	float vx, vy, vz;
	vx = 2.0f * (q_acc[1] * q_acc[3] - q_acc[0] * q_acc[2]);
	vy = 2.0f * (q_acc[0] * q_acc[1] + q_acc[2] * q_acc[3]);
	vz = 2.0f * (0.5f - q_acc[1]*q_acc[1] - q_acc[2]*q_acc[2]);

	float M[6] = { vx, vy, vz, Magnetometer[0], Magnetometer[1], Magnetometer[2] };

	if (fabs(1.0f/norm_m - std_norm_m) < threshold_m)
	{
		acc_mag(M, q_acc, q_gravity_mag);

		q_est[0] = (1.0f - gain_m)*q_acc[0] + gain_m*q_gravity_mag[0];
		q_est[1] = (1.0f - gain_m)*q_acc[1] + gain_m*q_gravity_mag[1];
		q_est[2] = (1.0f - gain_m)*q_acc[2] + gain_m*q_gravity_mag[2];
		q_est[3] = (1.0f - gain_m)*q_acc[3] + gain_m*q_gravity_mag[3];

		norm_q = invSqrt(q_est[0] * q_est[0] + q_est[1] * q_est[1] + q_est[2] * q_est[2] + q_est[3] * q_est[3]);
		q_est[0] *= norm_q;
		q_est[1] *= norm_q;
		q_est[2] *= norm_q;
		q_est[3] *= norm_q;
	}
	else
	{
		q_est[0] = q_acc[0];
		q_est[1] = q_acc[1];
		q_est[2] = q_acc[2];
		q_est[3] = q_acc[3];
	}

	
	



	return;

}


int main()
{

	float acc[3], gyro[3], mag[3], euler_true[3],time,dt,last_q[4],q_est[4];

	q_est[0] = 1.0f;
	q_est[1] = 0.0f;
	q_est[2] = 0.0f;
	q_est[3] = 0.0f;

	last_q[0] = 1.0f;
	last_q[1] = 0.0f;
	last_q[2] = 0.0f;
	last_q[3] = 0.0f;

	dt = 1.0f / 500.0f;

	long length = 0;
	std::string in;

	while (std::getline(input1, in))
		++length;

	for (int i = 0; i < length; ++i)
	{
		input2		   >> euler_true[0] >> euler_true[1] >> euler_true[2]
					   >> acc[0] >> acc[1] >> acc[2]
					   >> gyro[0] >> gyro[1] >> gyro[2]
			           >> mag[0] >> mag[1] >> mag[2] >>dt;

		acc[0] *= -1.0f;
		acc[2] *= -1.0f;

		gyro[0] *= -1.0f;
		gyro[2] *= -1.0f;

		mag[0] *= -1.0f;
		mag[2] *= -1.0f;


		FCF_gyro_acc_mag(gyro, acc, mag, last_q, 0.005f, 0.01f, 0.285f, 0.5f, dt, q_est);

		last_q[0] = q_est[0];
		last_q[1] = q_est[1];
		last_q[2] = q_est[2];
		last_q[3] = q_est[3];

		output << q_est[0]<<" "<< q_est[1]<<" " << q_est[2]<<" " << q_est[3] << std::endl;



	}

	return(0);
}



