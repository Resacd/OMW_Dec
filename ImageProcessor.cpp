#include "ImageProcessor.h"
#include <iostream>
/*
#include <opencv2/highgui.hpp>

ImageProcessor::ImageProcessor(const std::string& imageFile) : imageFile(imageFile) {}

void ImageProcessor::displayImage() {
    cv::Mat originalImage = cv::imread(imageFile);

    if (originalImage.empty()) {
        std::cerr << "Error: Could not open image." << std::endl;
        return;
    }

    cv::imshow("Original Image", originalImage);

    cv::waitKey(0);
}


void ImageProcessor::showContinuousBlackLine() {
    cv::Mat originalImage = cv::imread(imageFile);

    if (originalImage.empty()) {
        std::cerr << "Error: Could not open image." << std::endl;
        return;
    }

    cv::Mat grayImage;
    cv::cvtColor(originalImage, grayImage, cv::COLOR_BGR2GRAY);

    cv::Mat thresholded;
    cv::threshold(grayImage, thresholded, 100, 255, cv::THRESH_BINARY);

    std::vector<std::vector<cv::Point>> contours;
    cv::findContours(thresholded, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);

    cv::Mat reshapedImage = originalImage.reshape(1, originalImage.total());
    reshapedImage.convertTo(reshapedImage, CV_32F);

    int k = 3;
    cv::Mat labels, centers;
    cv::kmeans(reshapedImage, k, labels, cv::TermCriteria(cv::TermCriteria::EPS + cv::TermCriteria::MAX_ITER, 100, 0.2), 3, cv::KMEANS_RANDOM_CENTERS, centers);

    labels = labels.reshape(0, originalImage.rows);

    cv::Mat resultImage(originalImage.size(), CV_8UC3, cv::Scalar(255, 255, 255));

    for (size_t i = 0; i < contours.size(); ++i) {
        double area = cv::contourArea(contours[i]);
        double perimeter = cv::arcLength(contours[i], true);
        double circularity = 4 * CV_PI * area / (perimeter * perimeter);
        cv::Scalar meanColor = cv::mean(originalImage, thresholded);
        int clusterLabel = static_cast<int>(labels.at<int>(i, 0));

        if (clusterLabel == 2 && circularity > 0.6) {
            cv::drawContours(resultImage, contours, static_cast<int>(i), cv::Scalar(0, 0, 0), 2);
        }
    }

    cv::imshow("Continuous Black Line", resultImage);

    cv::waitKey(0);
}
*/
