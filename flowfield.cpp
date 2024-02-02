/*
#include <SFML/Graphics.hpp>
#include <vector>
#include <iostream>

const double D = 0.1;  // Molecular diffusion coefficient
const double u = 1.0;  // Fluid velocity
const double dx = 0.1; // Spatial step size
const double dt = 0.01; // Time step size
const double tMax = 2.0; // Maximum simulation time
const double initialConcentration = 1.0;

using namespace sf;

int main() {
    const int numX = static_cast<int>(1.0 / dx) + 1;
    const int numT = static_cast<int>(tMax / dt);

    std::vector<std::vector<double>> concentration(numT, std::vector<double>(numX, 0.0));

    concentration[0][numX / 4] = initialConcentration;

    sf::RenderWindow window(sf::VideoMode(), "Solute Dispersion");

    // Create a grid of vector fields
    sf::VertexArray vectors;
    //vectors.setPrimitiveType(sf::Lines);
    vectors.resize(numX * numT * 2);

    // Perform time-stepping to simulate dispersion
    for (int t = 0; t < numT - 1; ++t) {
        for (int x = 1; x < numX - 1; ++x) {
            concentration[t + 1][x] = concentration[t][x] + D * (concentration[t][x + 1] - 2 * concentration[t][x] + concentration[t][x - 1]) / (dx * dx) - u * (concentration[t][x + 1] - concentration[t][x - 1]) / (2 * dx) * dt;

            // Calculate the vectors for the vector field
            double concentrationChange = concentration[t + 1][x] - concentration[t][x];
            double xVelocity = concentrationChange * u / dx;
            
            vectors[2 * (t * numX + x)].position = sf::Vector2f(x * 10, t * 10);
            vectors[2 * (t * numX + x) + 1].position = sf::Vector2f(x * 10 + xVelocity * 10, t * 10);
        }
    }

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        window.clear(sf::Color::White);

        // Draw the vector field
        window.draw(vectors);

        window.display();
    }

    return 0;
}
*/

#include <iostream>
#include <vector>

const double D = 0.1;  // Molecular diffusion coefficient
const double u = 1.0;  // Fluid velocity
const double dx = 0.1;  // Spatial step size
const double dt = 0.01; // Time step size
const double tMax = 2.0; // Maximum simulation time
const double initialConcentration = 1.0;

int main() {
    // Calculate the number of spatial steps
    const int numX = static_cast<int>(1.0 / dx) + 1;

    // Calculate the number of time steps
    const int numT = static_cast<int>(tMax / dt);

    // Initialize concentration array
    std::vector<std::vector<double> > concentration(numT, std::vector<double>(numX, 0.0));

    // Set initial concentration at a specific location
    concentration[0][numX / 2] = initialConcentration;

    // Perform time-stepping to simulate dispersion
    for (int t = 0; t < numT - 1; ++t) {
        for (int x = 1; x < numX - 1; ++x) {
            // Calculate concentration change using the advection-dispersion equation
            concentration[t + 1][x] = concentration[t][x] + D * (concentration[t][x + 1] - 2 * concentration[t][x] + concentration[t][x - 1]) / (dx * dx) - u * (concentration[t][x + 1] - concentration[t][x - 1]) / (2 * dx) * dt;
        }
    }

    // Print the concentration profile at different time steps
    for (int t = 0; t < numT; t += numT / 10) {
        std::cout << "Time t = " << t * dt << ":\n";
        for (int x = 0; x < numX; ++x) {
            std::cout << abs(concentration[t][x]) << " ";
        }
        std::cout << "\n";
    }

    return 0;
}
