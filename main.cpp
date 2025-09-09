#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>

struct Particle {
    double mass;
    double x, y, z;
    double vx, vy, vz;
    double fx, fy, fz;
};

const double G = 6.67430e-11; // Gravitational constant

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cout << "Usage: " << argv[0] << " input.tsv dt steps dumpInterval\n";
        return 1;
    }

    std::string inputFile = argv[1];
    double dt = std::stod(argv[2]);
    int steps = std::stoi(argv[3]);
    int dumpInterval = std::stoi(argv[4]);

    std::ifstream file(inputFile);
    if (!file) {
        std::cerr << "Error: cannot open file " << inputFile << "\n";
        return 1;
    }

    // Read all numbers from the file
    std::vector<double> values;
    double val;
    while (file >> val) {
        values.push_back(val);
    }

    if (values.size() % 7 != 0) {
        std::cerr << "Error: input file does not have 7 numbers per particle.\n";
        return 1;
    }

    int n = values.size() / 7;
    std::vector<Particle> particles(n);

    // Assign numbers to particles
    for (int i = 0; i < n; ++i) {
        particles[i].mass = values[i * 7 + 0];
        particles[i].x    = values[i * 7 + 1];
        particles[i].y    = values[i * 7 + 2];
        particles[i].z    = values[i * 7 + 3];
        particles[i].vx   = values[i * 7 + 4];
        particles[i].vy   = values[i * 7 + 5];
        particles[i].vz   = values[i * 7 + 6];
    }

    // Main simulation loop
    for (int step = 0; step < steps; ++step) {
        for (auto &p : particles) p.fx = p.fy = p.fz = 0;

        // Compute pairwise gravitational forces
        for (int i = 0; i < n; ++i) {
            for (int j = i+1; j < n; ++j) {
                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;
                double dz = particles[j].z - particles[i].z;
                double dist2 = dx*dx + dy*dy + dz*dz + 1e-10;
                double dist = std::sqrt(dist2);
                double F = G * particles[i].mass * particles[j].mass / dist2;

                double Fx = F * dx / dist;
                double Fy = F * dy / dist;
                double Fz = F * dz / dist;

                particles[i].fx += Fx;
                particles[i].fy += Fy;
                particles[i].fz += Fz;

                particles[j].fx -= Fx;
                particles[j].fy -= Fy;
                particles[j].fz -= Fz;
            }
        }

        // Update velocities and positions
        for (auto &p : particles) {
            p.vx += p.fx / p.mass * dt;
            p.vy += p.fy / p.mass * dt;
            p.vz += p.fz / p.mass * dt;

            p.x += p.vx * dt;
            p.y += p.vy * dt;
            p.z += p.vz * dt;
        }

        // Output at dump intervals
        if (step % dumpInterval == 0) {
            std::ofstream out("output.tsv", std::ios::app); // append
            out << n;
            for (auto &p : particles) {
                out << "\t" << p.mass
                    << "\t" << p.x << "\t" << p.y << "\t" << p.z
                    << "\t" << p.vx << "\t" << p.vy << "\t" << p.vz
                    << "\t" << p.fx << "\t" << p.fy << "\t" << p.fz;
            }
            out << "\n";
        }
    }

    return 0;
}