#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

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
    int n;
    file >> n;
    std::vector<Particle> particles(n);

    for (int i = 0; i < n; i++) {
        file >> particles[i].mass
             >> particles[i].x >> particles[i].y >> particles[i].z
             >> particles[i].vx >> particles[i].vy >> particles[i].vz;
    }

    // Main simulation loop
    for (int step = 0; step < steps; step++) {
        // Reset forces
        for (auto &p : particles) p.fx = p.fy = p.fz = 0;

        // Compute pairwise gravitational forces
        for (int i = 0; i < n; i++) {
            for (int j = i+1; j < n; j++) {
                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;
                double dz = particles[j].z - particles[i].z;
                double dist2 = dx*dx + dy*dy + dz*dz + 1e-10; // avoid divide by 0
                double dist = std::sqrt(dist2);
                double F = G * particles[i].mass * particles[j].mass / dist2;

                double Fx = F * dx / dist;
                double Fy = F * dy / dist;
                double Fz = F * dz / dist;

                particles[i].fx += Fx;
                particles[i].fy += Fy;
                particles[i].fz += Fz;

                particles[j].fx -= Fx; // Newton's third law
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
            std::ofstream out("output.tsv");
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