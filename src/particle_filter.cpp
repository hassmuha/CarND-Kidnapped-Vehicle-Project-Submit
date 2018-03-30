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

	num_particles = 100; //TODO need to be changed

	Particle Particle_t;

	// This line creates a normal (Gaussian) distribution for x, y theta
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	//default random generator for each distribution
	default_random_engine gen_x;
	default_random_engine gen_y;
	default_random_engine gen_theta;

	for (unsigned int i = 0; i < num_particles; i++){
		Particle_t.id = i;
		Particle_t.x = dist_x(gen_x);
		Particle_t.y = dist_y(gen_x);
		Particle_t.theta = dist_theta(gen_x);
		Particle_t.weight = 1.0;
		particles.push_back(Particle_t);
		weights.push_back(1.0);
		cout<< "Particle_t.x " << Particle_t.x << endl;
		cout<< "Particle_t.y " << Particle_t.y << endl;
		cout<< "Particle_t.theta " << Particle_t.theta << endl;
	}
	cout<< "Total No of intialized samples" << particles.size() << endl;
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/


	//default random generator for each distribution
	default_random_engine gen_x;
	default_random_engine gen_y;
	default_random_engine gen_theta;
	//temporary values
	double x_pred;
  double y_pred;
  double theta_pred;


	for (unsigned int i = 0; i < num_particles; i++){

		// if yaw rate not equal to zero
		if (abs(yaw_rate) > 1e-5){
			x_pred = particles[i].x + (velocity/yaw_rate) * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			y_pred = particles[i].y + (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			theta_pred = particles[i].theta + yaw_rate*delta_t;
		} else {
			x_pred = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			y_pred = particles[i].y + velocity*delta_t*sin(particles[i].theta);
			theta_pred = particles[i].theta;
		}

		// This line creates a normal (Gaussian) distribution for x, y theta
		normal_distribution<double> dist_x(x_pred, std_pos[0]);
		normal_distribution<double> dist_y(y_pred, std_pos[1]);
		normal_distribution<double> dist_theta(theta_pred, std_pos[2]);

		particles[i].x = dist_x(gen_x);
		particles[i].y = dist_y(gen_x);
		particles[i].theta = dist_theta(gen_x);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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

	double sig_x = std_landmark[0];
	double sig_y = std_landmark[1];
	//calculate normalization term
	double gauss_norm= (1/(2 * M_PI * sig_x * sig_y));

	double x_p, y_p, theta, x_lm,y_lm , x_c, y_c ,x_m, y_m;
	int id_i;
	double mu_x, mu_y, exponent, P_m;

	//cout << "Total No of landmarks = " << map_landmarks.landmark_list.size() << endl;
	for (unsigned int i = 0; i < num_particles; i++){
		//Transforming the Vehicle Coordinate system observation to map coordinate system
		x_p = particles[i].x; //x position of particle in map coordinate
		y_p = particles[i].y; //y position of particle in map coordinate
		theta = particles[i].theta;

		// loop to find landmarks which can be covered by the particle current position and it's Range
		std::vector<int> landmarks_ids;
		for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++){
			x_lm = map_landmarks.landmark_list[j].x_f;
			y_lm = map_landmarks.landmark_list[j].y_f;
			//calcucate the euclidian distance and check with the range
			if (sqrt((x_lm-x_p)*(x_lm-x_p) + (y_lm-y_p)*(y_lm-y_p)) < sensor_range){
			//if (1){ //TODO Change this afterwards
				landmarks_ids.push_back(j);
			}
		}
		//cout << "No of landmarks = " << landmarks_ids.size() << " vs No of observations = " << observations.size()  << endl;


		std::vector<int> associations;
		std::vector<double> sense_x;
		std::vector<double> sense_y;

		// loop over observations
		//cout << "No of observations = "<< observations.size() << endl;
		for (unsigned int j = 0; j < observations.size(); j++){
			x_c = observations[j].x;
			y_c = observations[j].y;


			// calculating the observations in transformed coordinates
			// transform to map x coordinate
			x_m= x_p + (cos(theta) * x_c) - (sin(theta) * y_c);
			// transform to map y coordinate
			y_m= y_p + (sin(theta) * x_c) + (cos(theta) * y_c);


			//Finding the landmark which the predicted observation corresponds to
			double min_dist = 1e100;
			int min_idx = 0;
			for (unsigned int k = 0; k < landmarks_ids.size(); k++){
				x_lm = map_landmarks.landmark_list[landmarks_ids[k]].x_f;
				y_lm = map_landmarks.landmark_list[landmarks_ids[k]].y_f;
				//calcucate the euclidian distance
				double min_dist_t = sqrt((x_lm-x_m)*(x_lm-x_m) + (y_lm-y_m)*(y_lm-y_m));
				if (min_dist_t < min_dist){
					min_idx = k;
					min_dist = min_dist_t;
				}
			}
			//predicted measurement of landmark
			mu_x = map_landmarks.landmark_list[landmarks_ids[min_idx]].x_f;
			mu_y = map_landmarks.landmark_list[landmarks_ids[min_idx]].y_f;
			id_i = map_landmarks.landmark_list[landmarks_ids[min_idx]].id_i;
			associations.push_back(id_i);
			sense_x.push_back(x_m);
			sense_y.push_back(y_m);


			//calculate the probability for calculating weights
			// calculate exponent
			exponent= (pow((x_m - mu_x),2))/(2 * pow(sig_x,2)) + (pow((y_m - mu_y),2))/(2 * pow(sig_y,2));
			// calculate weight using normalization terms and exponent
			P_m= gauss_norm * exp(-exponent);
			weights[i] = weights[i]*P_m;
		}
		particles[i] = SetAssociations(particles[i],associations,sense_x,sense_y);
		particles[i].weight = weights[i];
	}

	int N = num_particles;
	//Find the weight sum for normalization and find the maxmimum weight
	double sum_w = 0.0;
	double max_w = 0.0;
	for (unsigned int i = 0; i < N; i++){
		sum_w += weights[i];
		if (max_w < weights[i]){
			max_w = weights[i];
		}
	}
	//cout << "Debug : Sum of Weight" << sum_w << endl;
	//cout << "Debug : Max Weight" << max_w << endl;

	//loop for normalization
	for (unsigned int i = 0; i < N; i++){
		weights[i] = weights[i]/sum_w;
		particles[i].weight = weights[i];
		//cout << "Debug : normalized weight " << weights[i] << endl;
	}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	int N = num_particles;
	default_random_engine gen;

	// some temporary vectors
	std::vector<Particle> particles_upd;
	std::vector<double> weights_upd;

	discrete_distribution<double> dist_rand (weights.begin(), weights.end());

	//loop for resampling
	for (unsigned int i = 0; i < N; i++){
		int index = dist_rand(gen);
		//cout << "Debug : Index Resamples " << index << endl;
    particles_upd.push_back(particles[index]);
		weights_upd.push_back(1.0); //new weight assign to 1.0
	}
	weights = weights_upd;
	particles = particles_upd;


}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
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
