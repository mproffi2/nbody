#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

const double G = 6.674e-11;   // gravitational constant
const double softening = 1e-9; // avoid divide-by-zero

// Particle structure
struct Particle {
    double mass;
    double x, y, z;
    double vx, vy, vz;
    double fx, fy, fz;
};

// Load from file (e.g. solar.tsv)
vector<Particle> loadFromFile(const string& filename) {
    ifstream file(filename);
    vector<Particle> particles;
    int n;
    if (file >> n) {
        particles.resize(n);
        for (int i = 0; i < n; i++) {
            Particle p;
            file >> p.mass >> p.x >> p.y >> p.z
                 >> p.vx >> p.vy >> p.vz
                 >> p.fx >> p.fy >> p.fz;
            particles[i] = p;
        }
    }
    return particles;
}

// Compute gravitational forces
void computeForces(vector<Particle>& particles) {
    for (auto& p : particles) {
        p.fx = p.fy = p.fz = 0.0;
    }
    int n = particles.size();
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double dx = particles[j].x - particles[i].x;
            double dy = particles[j].y - particles[i].y;
            double dz = particles[j].z - particles[i].z;
            double distSqr = dx*dx + dy*dy + dz*dz + softening;
            double dist = sqrt(distSqr);

            double F = (G * particles[i].mass * particles[j].mass) / distSqr;
            double fx = F * dx / dist;
            double fy = F * dy / dist;
            double fz = F * dz / dist;

            particles[i].fx += fx;
            particles[i].fy += fy;
            particles[i].fz += fz;

            particles[j].fx -= fx;
            particles[j].fy -= fy;
            particles[j].fz -= fz;
        }
    }
}

// Update positions and velocities
void update(vector<Particle>& particles, double dt) {
    for (auto& p : particles) {
        double ax = p.fx / p.mass;
        double ay = p.fy / p.mass;
        double az = p.fz / p.mass;

        p.vx += ax * dt;
        p.vy += ay * dt;
        p.vz += az * dt;

        p.x += p.vx * dt;
        p.y += p.vy * dt;
        p.z += p.vz * dt;
    }
}

// Output in TSV format
void printState(const vector<Particle>& particles) {
    cout << particles.size();
    for (auto& p : particles) {
        cout << "\t" << p.mass << "\t" << p.x << "\t" << p.y << "\t" << p.z
             << "\t" << p.vx << "\t" << p.vy << "\t" << p.vz
             << "\t" << p.fx << "\t" << p.fy << "\t" << p.fz;
    }
    cout << endl;
}

// Main loop
int main(int argc, char* argv[]) {
    if (argc < 5) {
        cerr << "Usage: " << argv[0] << " input.tsv dt steps dumpInterval\n";
        return 1;
    }

    string filename = argv[1];
    double dt = stod(argv[2]);
    int steps = stoi(argv[3]);
    int dumpInterval = stoi(argv[4]);

    vector<Particle> particles = loadFromFile(filename);

    for (int step = 0; step < steps; step++) {
        computeForces(particles);
        update(particles, dt);
        if (step % dumpInterval == 0) {
            printState(particles);
        }
    }
    return 0;
}