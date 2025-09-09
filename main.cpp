#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>

// Structure to hold particle data
struct Particle {
    double mass = 0;
    double x = 0, y = 0, z = 0;
    double vx = 0, vy = 0, vz = 0;
    double fx = 0, fy = 0, fz = 0;
};

// Gravitational constant
const double G = 6.67430e-11;

// Function to safely read a number, returns 0 if missing
double safeRead(std::istringstream &ss) {
    double val = 0;
    ss >> val;
    if (ss.fail()) {
        ss.clear();
        val = 0;
    }
    return val;
}

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
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << inputFile << "\n";
        return 1;
    }

    int n;
    file >> n;
    std::vector<Particle> particles(n);

    // Read particle data safely
    for (int i = 0; i < n; i++) {
        std::string line;
        std::getline(file, line); // read rest of the line
        if(line.empty()) std::getline(file, line);
        std::istringstream ss(line);
        particles[i].mass = safeRead(ss);
        particles[i].x    = safeRead(ss);
        particles[i].y    = safeRead(ss);
        particles[i].z    = safeRead(ss);
        particles[i].vx   = safeRead(ss);
        particles[i].vy   = safeRead(ss);
        particles[i].vz   = safeRead(ss);
    }

    // Main simulation loop
    for (int step = 0; step < steps; step++) {
        // Reset forces
        for (auto &p : particles) p.fx = p.fy = p.fz = 0;

        // Compute gravitational forces
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

                particles[j].fx -= Fx;
                particles[j].fy -= Fy;
                particles[j].fz -= Fz;
            }
        }

        // Update velocities and positions safely
        for (auto &p : particles) {
            if (p.mass != 0) {
                p.vx += p.fx / p.mass * dt;
                p.vy += p.fy / p.mass * dt;
                p.vz += p.fz / p.mass * dt;
            }
            p.x += p.vx * dt;
            p.y += p.vy * dt;
            p.z += p.vz * dt;
        }

        // Dump output at intervals
        if (step % dumpInterval == 0) {
            std::ofstream out("output.tsv", step == 0 ? std::ios::trunc : std::ios::app);
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