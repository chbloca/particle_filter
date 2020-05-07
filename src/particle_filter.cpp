/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles

  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  std::default_random_engine gen;

  weights.resize(num_particles, 1.0);

  for(int i = 0; i < num_particles; i++){
      Particle particle;
      particle.id = i;
      particle.x = dist_x(gen);
      particle.y = dist_y(gen);
      particle.theta = dist_theta(gen);
      particle.weight = 1.0;
      particles.push_back(particle);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
    double pred_x;
    double pred_y;
    double pred_theta;
    std::default_random_engine gen;

    for(int i = 0; i < num_particles; i++){
        if(fabs(yaw_rate < 0.0001)){
            pred_x = particles[i].x + velocity * cos(particles[i].theta) * delta_t;
            pred_y = particles[i].y + velocity * sin(particles[i].theta) * delta_t;
            pred_theta = particles[i].theta;
        }else{
            pred_x = particles[i].x + velocity / yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
            pred_y = particles[i].y + velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta - yaw_rate * delta_t));
            pred_theta = particles[i].theta + yaw_rate * delta_t;
        }

        normal_distribution<double> dist_x(pred_x, std_pos[0]);
        normal_distribution<double> dist_y(pred_y, std_pos[1]);
        normal_distribution<double> dist_theta(pred_theta, std_pos[2]);

        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
    }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
    double min_distance = std::numeric_limits<double>::max();
    int nearest_id;
    for(int i = 0; i < observations.size(); i++){
        for(int j = 0; j < predicted.size(); j++){
            double distance = dist(observations[i].x, predicted[j].x, observations[i].y, predicted[j].y);
            if(distance < min_distance){
                min_distance = distance;
                nearest_id = predicted[j].id;
            }
        }
        observations[i].id = nearest_id;
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
    double weight_total = 0.0;
    for(int i = 0; i < num_particles; i++){

        // Apply homogeneous transformation to conver from car coordinate system to map coordinate system
        std::vector<LandmarkObs> transformed_observations;
        for(int j = 0; j < observations.size(); j++){
            LandmarkObs tmp;
            tmp.x = cos(particles[i].theta) * observations[j].x - sin(particles[i].theta) * observations[j].x;
            tmp.y = sin(particles[i].theta) * observations[j].y + cos(particles[i].theta) * observations[j].y;
            transformed_observations.push_back(tmp);
        }
        // Discard landmarks out of sensor range
        std::vector<LandmarkObs> landmarks_in_range;
        for(int k = 0; k < map_landmarks.landmark_list.size(); k++){
            double range = dist(particles[i].x, particles[i].y,
                                map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);
            if(range < sensor_range){
                LandmarkObs tmp;
                tmp.x = map_landmarks.landmark_list[k].x_f;
                tmp.y = map_landmarks.landmark_list[k].y_f;
                tmp.id = map_landmarks.landmark_list[k].id_i;
                landmarks_in_range.push_back(tmp);
            }
        }

        dataAssociation(landmarks_in_range, transformed_observations);
        double sig_x = std_landmark[0];
        double sig_y = std_landmark[1];
        double exponent;

        double gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);
        for(int l = 0; transformed_observations.size(); l++){
            double x = transformed_observations[l].x;
            double y = transformed_observations[l].y;

            for(int m = 0; map_landmarks.landmark_list.size(); m++){
                if(map_landmarks.landmark_list[m].id_i == transformed_observations[l].id){
                    double mu_x = landmarks_in_range[m].x;
                    double mu_y = landmarks_in_range[m].y;

                    exponent = pow(x - mu_x, 2) / (2 * sig_x * sig_x) +
                                      pow(y - mu_y, 2) / (2 * sig_y * sig_y);

                    particles[i].weight *= gauss_norm * exp(- exponent);
                }
            }
            weights[i] = gauss_norm * exp(- exponent);
        }
        weight_total += particles[i].weight;
    }

    for(int i = 0; i < particles.size(); i++){
        particles[i].weight /= weight_total;
        weights[i] = particles[i].weight;
    }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */


}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
