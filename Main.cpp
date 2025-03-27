#include <iostream>
#include <armadillo>
#include <cmath>

using namespace arma;
using namespace std;

int main() {
    // Parameters
    int N = 1000;            // Number of grid points
    double x_min = -10.0;    // Left boundary (units of sqrt(hbar/mω))
    double x_max = 10.0;     // Right boundary
    double dx = (x_max - x_min) / (N - 1);  // Spatial step
    vec x = linspace(x_min, x_max, N);      // Position grid

    // Hamiltonian matrix (finite difference method)
    mat H = zeros<mat>(N, N);
    double hbar = 1.0;       // Natural units (hbar = m = ω = 1)
    double m = 1.0;
    double omega = 1.0;

    // Kinetic energy (T = -ħ²/2m d²/dx²)
    double T = -hbar * hbar / (2 * m * dx * dx);
    // Potential energy (V = 0.5 m ω² x²)
    vec V = 0.5 * m * omega * omega * x % x;  // Element-wise square

    // Construct Hamiltonian
    for (int i = 1; i < N - 1; i++) {
        H(i, i - 1) = T;
        H(i, i) = -2 * T + V(i);  // Diagonal: kinetic + potential
        H(i, i + 1) = T;
    }

    // Boundary conditions (ψ → 0 at x = ±∞)
    H(0, 0) = 1e30;       // Approximate ∞ potential
    H(N - 1, N - 1) = 1e30;   // Approximate ∞ potential

    // Solve eigenvalue problem (Hψ = Eψ)
    vec energies;
    mat wavefunctions;
    eig_sym(energies, wavefunctions, H);

    // Print first 5 energy levels (should be E_n = n + 0.5 in natural units)
    cout << "Quantum Harmonic Oscillator Energy Levels:" << endl;
    for (int n = 0; n < 5; n++) {
        cout << "E_" << n << " = " << energies(n) << " (Expected: " << n + 0.5 << ")" << endl;
    }

    // Save wavefunctions to file
    wavefunctions.save("qho_wavefunctions.txt", raw_ascii);

    return 0;
}