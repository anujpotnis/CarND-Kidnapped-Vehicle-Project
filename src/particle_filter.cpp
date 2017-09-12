/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    
    std::default_random_engine gen;
    
    
    
    num_particles = 15;
    
    
    std::normal_distribution<double> dist_x(x, std[0]);
    std::normal_distribution<double> dist_y(y, std[1]);
    std::normal_distribution<double> dist_theta(theta, std[2]);
    
    for (int i=0; i < num_particles; i++) {
        Particle particle;
        particle.id = i;
        particle.x = dist_x(gen);
        particle.y = dist_y(gen);
        particle.theta = dist_theta(gen);
        particle.weight = 1;
        
        particles.push_back(particle);
        weights.push_back(1);
    }
    
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/
    
    
    default_random_engine gen;
    for (int i = 0; i < num_particles; i++) {
        double new_x;
        double new_y;
        double new_theta;
        
        if (abs(yaw_rate) < 0.001) {
            new_x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
            new_y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
            new_theta = particles[i].theta;
        }
        else {
            new_x = particles[i].x + velocity/yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
            new_y = particles[i].y + velocity/yaw_rate * (-cos(particles[i].theta + yaw_rate*delta_t) + cos(particles[i].theta));
            new_theta = particles[i].theta + yaw_rate*delta_t;
        }
        
        std::normal_distribution<double> dist_x(new_x, std_pos[0]);
        std::normal_distribution<double> dist_y(new_y, std_pos[1]);
        std::normal_distribution<double> dist_theta(new_theta, std_pos[2]);
        
        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.
    
    for (int i = 0; i < observations.size(); i++) {
        int map_id = -1;
        double min_dist = numeric_limits<double>::max();
        
        for (int j = 0; j < predicted.size(); j++) {
            double cur_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
            
            if (cur_dist < min_dist) {
                min_dist = cur_dist;
                map_id = predicted[j].id;
            }
        }
        observations[i].id = map_id;
        
    }
}



void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   std::vector<LandmarkObs> observations, Map map_landmarks) {
    // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    //   according to the MAP'S coordinate system. You will need to transform between the two systems.
    //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
    //   The following is a good resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    //   and the following is a good resource for the actual equation to implement (look at equation
    //   3.33
    //   http://planning.cs.uiuc.edu/node99.html
    
    
    for (int i = 0; i < num_particles; i++) {
        vector<LandmarkObs> inRange_observations;
        
        for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
            // assigning indexed values of vector<> speeds up operations !!
            double lm_x = map_landmarks.landmark_list[j].x_f;
            double lm_y = map_landmarks.landmark_list[j].y_f;
            int lm_id = map_landmarks.landmark_list[j].id_i;
            // consider only landmarks within sensor range
            // DO NOT USE dist(). It is SLOW
            //if (dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[i].x_f, map_landmarks.landmark_list[i].y_f) <= sensor_range) {
            if (fabs(particles[i].x - lm_x) <= sensor_range && fabs(particles[i].y - lm_y) <= sensor_range) {
                inRange_observations.push_back(LandmarkObs{lm_id, lm_x, lm_y});
            }
        }
        
        // take observations of one particle, convert from vehicle to map co-ordinates
        vector<LandmarkObs> t_observations;
        double p_x = particles[i].x;
        double p_y = particles[i].y;
        double p_theta = particles[i].theta;
        
        for (int j = 0; j<observations.size(); j++) {
            double t_x = cos(p_theta) * observations[j].x - sin(p_theta) * observations[j].y + p_x;
            double t_y = sin(p_theta) * observations[j].x + cos(p_theta) * observations[j].y + p_y;
            t_observations.push_back(LandmarkObs{observations[j].id, t_x, t_y});
        }
        
        // associate the transformed observation to one of the in range observations (nearest)
        //std::cout << inRange_observations.size() << " " << t_observations.size() << std::endl;
        dataAssociation(inRange_observations, t_observations);
        //std::cout << inRange_observations.size() << " " << t_observations.size() << std::endl;
        
        // Reinitialize weights after particles have been updated
        particles[i].weight = 1.0;
        
        double s_x = std_landmark[0];
        double s_y = std_landmark[1];
        
        // pre-calculate denominator for speeding up computation
        const double gauss_d = (2*M_PI*s_x*s_y);
        const double gauss_d_x = 0.5 * s_x * s_x;
        const double gauss_d_y = 0.5 * s_x * s_x;
        
        for (int j = 0; j < t_observations.size(); j++) {
            double t_o_x, t_o_y, iRo_x, iRo_y;
            t_o_x = t_observations[j].x;
            t_o_y = t_observations[j].y;
            int associated_inRange_observation = t_observations[j].id;
            for (int k = 0; k < inRange_observations.size(); k++) {
                if (inRange_observations[k].id == associated_inRange_observation) {
                    iRo_x = inRange_observations[k].x;
                    iRo_y = inRange_observations[k].y;
                }
            }
            
            // Multivariate Gaussian
            double w_observation = (1/gauss_d) * exp(-(pow(iRo_x-t_o_x,2) / gauss_d_x + (pow(iRo_y-t_o_y,2) / gauss_d_y)));
            particles[i].weight *= w_observation;
        }
    }
}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    std::default_random_engine gen;
    // get all of the current weights
    // IMPORTANT (from Slack) : "You need to clear the `weights` vector before you start pushing the weights."
    vector<double> weights;
    for (int i = 0; i < num_particles; i++) {
        weights.push_back(particles[i].weight);
    }
    std::discrete_distribution<int> distribution(weights.begin(), weights.end());
    vector<Particle> resample_particles(num_particles);
    for (int i=0; i<num_particles; i++) {
        resample_particles[i] = (particles[distribution(gen)]);
    }
    particles = std::move(resample_particles);
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
    
    //Clear the previous associations
    particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();
    
    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
    
    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
    vector<int> v = best.associations;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
    vector<double> v = best.sense_x;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
    vector<double> v = best.sense_y;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
