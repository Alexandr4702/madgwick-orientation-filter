/*
 * Madgwick.h
 *
 *  Created on: Mar 26, 2021
 *      Author: alexandr
 */

#ifndef MADGWICK_H_
#define MADGWICK_H_

#include "Eigen/Eigen"

using Vec3 = Eigen::Vector3d;
using Quat = Eigen::Quaternion<double>;

/*
 *@brief: functions to work with quaternions
 *
*/
//----------------------------------------------------------------------------------------------------------------
Quat operator* (Quat& q1, double Scalar );
Quat operator* (Quat&& q1, double Scalar );

Quat operator* (double Scalar, Quat& q1 );
Quat operator* (double Scalar, Quat&& q1 );
//----------------------------------------------------------------------------------------------------------------
Quat operator/ (Quat& q1, double Scalar );
Quat operator/ (Quat&& q1, double Scalar );
Quat operator/ (double Scalar, Quat& q1 );
Quat operator/ (double Scalar, Quat&& q1 );
//----------------------------------------------------------------------------------------------------------------
Quat operator+ (Quat& q1, double Scalar );
Quat operator+ (Quat&& q1, double Scalar );
Quat operator+ (double Scalar, Quat& q1 );
Quat operator+ (double Scalar, Quat&& q1 );

Quat operator+ (Quat& q1, Quat& q2 );
Quat operator+ (Quat&& q1, Quat& q2 );
Quat operator+ (Quat& q1, Quat&& q2 );
Quat operator+ (Quat&& q1, Quat&& q2 );
//----------------------------------------------------------------------------------------------------------------
Quat operator- (Quat& q1, double Scalar );
Quat operator- (Quat&& q1, double Scalar );
Quat operator- (double Scalar, Quat& q1 );
Quat operator- (double Scalar, Quat&& q1 );

Quat operator- (Quat& q1, Quat& q2 );
Quat operator- (Quat&& q1, Quat& q2 );
Quat operator- (Quat& q1, Quat&& q2 );
Quat operator- (Quat&& q1, Quat&& q2 );
//----------------------------------------------------------------------------------------------------------------


class Madgwick_filter
{
public:
    /*!
     * \brief update all vectors have to be normalized otherwise it won't work correct
     * \param sun_ref
     * \param magn_ref
     * \param sun_meas
     * \param magn_meas
     * \param omega
     */
    void update(Vec3 sun_ref, Vec3 magn_ref, Vec3 sun_meas, Vec3 magn_meas, Vec3 omega, double dt);
    void update(Vec3 sun_ref, Vec3 sun_meas, Vec3 omega, double dt);


    Quat get_orientation();

    Vec3 getOmega_bias() const;
    void setOmega_bias(const Vec3 &value);
    /*!
     * \brief get_omega_from_quat
     * \param q
     * \param q_1
     * \param dt
     * \return
     */
    static Vec3 get_omega_from_quat(Quat q, Quat q_1, double dt);

private:
    /*!
     * \brief omega_integration
     * \param q_1
     * \param omega
     * \param dt
     * \return
     */
    Quat omega_integration (Quat q_1,Vec3 omega, double dt);

    /*!
     * \brief func
     * \param q quat
     * \param d ref
     * \param s sensors measurment
     * \return q^{-1} * d * q - s
     */
    Quat func(Quat q, Vec3 d, Vec3 s);

    Eigen::Matrix<double, 3, 4> Jacobian(Quat q, Vec3 d);

    /*!
     * \brief grad_function q^{-1} * d * q - s
     * \param q quat
     * \param d reference vector
     * \param s sensor measurments
     * \return J^T * f
     */
    Quat grad_function(Vec3 ref, Vec3 meas, Quat q);

    /*!
     * \brief returns q dervative from current quaternion and current omega
     * \param q_1
     * \param omega
     * \return
     */
    Quat q_dervative_omega (Quat q_1, Vec3 omega);

private:
	double betta = 0.1;
	Quat q_1 = Quat(1,0,0,0); //we should begin from quat from triad
	Vec3 omega_bias = Vec3(0,0,0);
};



#endif /* MADGWICK_H_ */
