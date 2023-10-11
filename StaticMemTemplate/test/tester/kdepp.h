#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "./kde.h"

namespace kdepp
{

    namespace kdemath
    {

        /// @brief Calculates covariance matrix.
        Eigen::MatrixXf covariance(std::vector<Eigen::VectorXf> const &data)
        {
            // TODO add checks that we have data at all?
            auto dim = data.front().size();
            Eigen::VectorXf mean = Eigen::VectorXf::Zero(dim);
            for (auto i : data)
            {
                mean += i;
            }
            mean = mean / (1.0 * data.size());

            Eigen::MatrixXf cov = Eigen::MatrixXf::Zero(dim, dim);
            for (auto x_j : data)
            {
                Eigen::MatrixXf tmp((x_j - mean) * (x_j - mean).transpose());
                cov += tmp;
            }
            cov = (1.0 / (data.size() - 1)) * cov; // should be n - 1?
            return cov;
        }

    } // namespace kdemath

    class Kdemv
    {
    public:
        Kdemv(std::vector<Eigen::VectorXf> data,
              std::string bandwidth_method = "scott")
            : data_(data), h_(Eigen::MatrixXf::Zero(data_.front().size(), data.front().size())), dim_(data_.front().size())
        {
            init_bandwidth(bandwidth_method);
            // debug h and h_inv_
            std::cout << "check h" << std::endl;
            std::cout << h_ << std::endl;
            std::cout << "check h_inv_" << std::endl;
            std::cout << h_inv_ << std::endl;
            pre_calculate_terms();
            std::cout << "pow_pi_term_ " << pow_pi_term_ << " h_pow_term_ " << h_pow_term_ << std::endl;
        };

        void set_bandwidth(Eigen::MatrixXf h)
        {
            h_ = h;
            h_inv_ = h_.inverse();
            pre_calculate_terms();
        };

        double eval(Eigen::VectorXf point)
        {
            double sum = 0;
            for (auto const &i : data_)
            {
                std::cout << "sample is" << std::endl;
                std::cout << i << std::endl;
                sum += kernel(point - i);
            }
            double n = 1.0 * data_.size();
            return sum / n;
        };

    private:
        using index_t = std::ptrdiff_t; // core guidelines

        double kernel(Eigen::VectorXf diff)
        {
            std::cout << "debug kernel diff " << std::endl;
            std::cout << diff << std::endl;
            double res = pow_pi_term_ * h_pow_term_ * std::exp((-1.0 / 2.0) * diff.transpose() * h_inv_ * diff);
            std::cout << "debug kernel res " << res << std::endl;
            return res;
        }

        void init_bandwidth(std::string bandwidth_method)
        {
            auto cov = kdemath::covariance(data_);
            std::cout << "cov matrix" << std::endl;
            std::cout << cov << std::endl;

            double n_term = std::pow(data_.size(), -1.0 / (dim_ + 4.0));
            std::cout << "n_term " << n_term << std::endl;

            if (bandwidth_method == "silverman")
            {
                double silverman_term = std::pow(4 / (dim_ + 2), -1.0 / (dim_ + 4.0));
                h_ = cov * (silverman_term * n_term * silverman_term * n_term);
            }
            else
            {
                // scott:
                h_ = cov * n_term * n_term;
            }
            h_inv_ = h_.inverse();
        }

        void pre_calculate_terms()
        {
            // for optimization do some pre calcs of power constants:
            pow_pi_term_ = std::pow(2 * kdemath::pi<double>(), -1.0 * dim_ / 2.0);
            h_pow_term_ = std::pow(h_.determinant(), -1.0 / 2.0);
        }

        std::vector<Eigen::VectorXf> data_;
        Eigen::MatrixXf h_;
        Eigen::MatrixXf h_inv_;

        // number of dimensions in problem:
        double dim_;

        // Pre calculated terms for efficiency only:
        double pow_pi_term_;
        double h_pow_term_;
    };

} // namespace kdepp