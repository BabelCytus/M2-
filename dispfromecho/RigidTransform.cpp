#include  <iostream>
#include <random>

#include <Eigen/Geometry>

#include "RigidTransform.h"

// Rotate a set of points aroud (center, axis) with an angle in degrees
void Rotation(const Points& points, Eigen::Vector3d center, Eigen::Vector3d axis, double angle, Points& newPoints)
{
    Eigen::Affine3d T = Eigen::Translation3d(center) * Eigen::AngleAxisd(angle * pi / 180, axis) 
                        * Eigen::Translation3d(-center);

    for (const auto& point : points) 
    {
        newPoints.push_back(Eigen::Vector3d{ T * point });
    }
}

// Translate a set of points
void Translation(const Points& points, const Eigen::Vector3d& t, Points& newPoints) 
{
    for (const auto& point : points) 
    {
        newPoints.push_back(point + t);
    }
}

// Estimation of rigid transform between two sets of points
// Arun's method ( https://jingnanshi.com/blog/arun_method_for_3d_reg.html )
RigidTransform Procrustes(const Points& pointsA, const Points& pointsB)
{
    Eigen::MatrixXd PA(pointsA.size(), 3), PB(pointsB.size(), 3);
    
    int i{};
    for (auto itA{ pointsA.cbegin() }, itB{ pointsB.cbegin() }; itA != pointsA.end(); ++itA, ++itB) {
        PA.row(i) = itA->transpose();
        PB.row(i) = itB->transpose();
        ++i;
    }

    // Pointset centers
    Eigen::Vector3d CA{ PA.colwise().mean() }, CB{ PB.colwise().mean() };
    //std::cout << CA << std::endl << CB << std::endl;

    // Pointset centering
    for (int i{}; i < PA.rows(); i++) PA.row(i) -= CA.transpose();
    for (int i{}; i < PB.rows(); i++) PB.row(i) -= CB.transpose();
    //std::cout << PA << std::endl << PB << std::endl;

    // Correlation matrix
    Eigen::Matrix3d K{ PB.transpose() * PA };
    //std::cout << "Correlation matrix:" << std::endl << K << std::endl;

    // Singular Value Decomposition of K (SVD)
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(K, Eigen::ComputeThinU | Eigen::ComputeThinV);
    //std::cout << svd.rank() << std::endl;
    //JacobiSVD<MatrixXd> svd(K, ComputeFullU | ComputeFullV);

    auto U{ svd.matrixU() };
    auto V{ svd.matrixV() };
    Eigen::DiagonalMatrix<double, 3> d(1, 1, V.determinant() * U.transpose().determinant());
    Eigen::Matrix3d RAB{ U * d * V.transpose() };
    //std::cout << "Rotation matrix:" << std::endl << RAB << std::endl;

    Eigen::Vector3d tAB{ CB - RAB * CA };
    //std::cout << "Translation vector:" << std::endl << tAB << std::endl;

    //for (int i{}; i < pointsA.size(); ++i)
    //    std::cout << pointsB[i].transpose() << "\t" << (RAB * pointsA[i] + tAB).transpose() << "\n";

    return { RAB, tAB };
}

std::tuple<RigidTransform, Indices, double> RANSAC(const Points& modelPoints, const Points& dataPoints, unsigned int sampleMinimalSize, unsigned int nIteration, double distanceThreshold, unsigned int nInliers)
{
    const double squaredDistanceThreshold{ distanceThreshold * distanceThreshold };
    Points dataSample;		    // Sample of n randomly chosen data points
    Points modelSample;         // Sample of n randomly chosen model points
    Indices inliers;            // Indices of inlier points

    Indices indices;	        // Indices of reshuffled points
    for (int i{}; i < dataPoints.size(); ++i) indices.push_back(i);

    // Rigid transform that best fits data (initialized to identity rigid transform
    RigidTransform bestModel{ Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero() };           
    Indices bestInliers;                // Inliers corresponding to the best model
    double bestError{ 10000.0 };        // Error between data and best model
    auto inliersNumber{ nInliers };

    std::random_device rd;
    std::mt19937 generator(rd());

    // Iterate up to nIteration
    for (unsigned int it{}; it < nIteration; ++it) 
    {
        std::shuffle(indices.begin(), indices.end(), generator);
        std::uniform_int_distribution<> dist(sampleMinimalSize, dataPoints.size());
        auto nSamplePoints{ dist(generator) };
        int n{};
        Points reshuffledModel, reshuffledData;
        for (auto ind{ indices.cbegin() }; ind != indices.cend(); ++ind) 
        {
            if (n < nSamplePoints) {
                dataSample.push_back(dataPoints[*ind]);
                modelSample.push_back(modelPoints[*ind]);
            }
            reshuffledData.push_back(dataPoints[*ind]);
            reshuffledModel.push_back(modelPoints[*ind]);
            ++n;
        }

        auto [R, t] { Procrustes(modelSample, dataSample) };
        //std::cout << "Model:\n" << R << std::endl << t << std::endl;

        // Apply rigid transform estimated on random sample all model points
        modelSample.clear();
        dataSample.clear();
        auto itM{ reshuffledModel.cbegin() }, itD{ reshuffledData.cbegin() };
        auto itI{ indices.cbegin() };
        while (itM < reshuffledModel.cend())
        {
            auto transformedPoint{ R * *itM + t };
            auto squaredDistance{ (transformedPoint - *itD).squaredNorm() };
            
            // Keep inliers in model and data samples
            if (squaredDistance < squaredDistanceThreshold)
            {
                modelSample.push_back(*itM);
                dataSample.push_back(*itD);
                inliers.push_back(*itI);
            }
            ++itM, ++itD, ++itI;
        }

        // Model is evaluated if number of inliers is higher than current minimum number of inliers
        if (dataSample.size() > inliersNumber) 
        {
            auto [R, t] { Procrustes(modelSample, dataSample) };

            double error{};
            for (auto itM{ modelSample.cbegin() }, itD{ dataSample.cbegin() }; itM < modelSample.cend(); ++itM, ++itD)
            {
                auto transformedPoint{ R * *itM + t };
                auto distance{ (transformedPoint - *itD).norm() };
                error += distance;
            }
            error /= modelSample.size();

            if (error < bestError) 
            {
                bestError = error;
                bestModel = { R, t };
                bestInliers = inliers;
                inliersNumber = modelSample.size();
            }
        }

        inliers.clear();
        modelSample.clear();
        dataSample.clear();
    }

    return { bestModel, bestInliers, bestError };
}
