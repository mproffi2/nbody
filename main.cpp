#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>  // for std::scientific

const double G = 6.674e-11;    // Gravitational constant
const double SOFTENING = 1e3;  // Softening factor to avoid singularities

// Particle structure
struct Particle {
    double mass;
    double x, y, z;
    double vx, vy, vz;
    double fx, fy, fz;
};

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " input_file dt steps dump_interval\n";
        return 1;
    }

    std::string filename = argv[1];
    double dt = std::stod(argv[2]);
    int steps = std::stoi(argv[3]);
    int dump_interval = std::stoi(argv[4]);

    // Open input file
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Cannot open file " << filename << "\n";
        return 1;
    }

    int N;
    infile >> N;  // number of particles
    std::vector<Particle> particles(N);

    // Read particle data from file
    for (int i = 0; i < N; ++i) {
        double mass, x, y, z, vx, vy, vz;
        infile >> mass >> x >> y >> z >> vx >> vy >> vz;
        particles[i] = {mass, x, y, z, vx, vy, vz, 0, 0, 0};

        // Skip the rest of the line (extra columns, if any)
        std::string dummy;
        std::getline(infile, dummy);
    }
    infile.close();

    // Open output file
    std::ofstream outfile("output.tsv");
    if (!outfile) {
        std::cerr << "Cannot open output.tsv\n";
        return 1;
    }

    outfile << std::scientific << std::setprecision(6); // readable scientific format

    // Main simulation loop
    for (int step = 0; step <= steps; ++step) {

        // Reset forces
        for (auto &p : particles)
            p.fx = p.fy = p.fz = 0.0;

        // Compute gravitational forces
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                double dx = particles[j].x - particles[i].x;
                double dy = particles[j].y - particles[i].y;
                double dz = particles[j].z - particles[i].z;
                double distSqr = dx*dx + dy*dy + dz*dz + SOFTENING*SOFTENING;
                double distSixth = distSqr * std::sqrt(distSqr);
                double F = G * particles[i].mass * particles[j].mass / distSixth;

                // Update forces on particle i
                particles[i].fx += F * dx;
                particles[i].fy += F * dy;
                particles[i].fz += F * dz;

                // Update forces on particle j (Newton's 3rd law)
                particles[j].fx -= F * dx;
                particles[j].fy -= F * dy;
                particles[j].fz -= F * dz;
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

        // Dump state to output file at intervals
        if (step % dump_interval == 0) {
            outfile << N;
            for (auto &p : particles) {
                outfile << "\t" << p.mass
                        << "\t" << p.x << "\t" << p.y << "\t" << p.z
                        << "\t" << p.vx << "\t" << p.vy << "\t" << p.vz
                        << "\t" << p.fx << "\t" << p.fy << "\t" << p.fz;
            }
            outfile << "\n";
        }
    }

    outfile.close();
    
    return 0;
}