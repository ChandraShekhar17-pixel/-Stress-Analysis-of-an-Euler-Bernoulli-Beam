# Euler-Bernoulli Beam Stress Analysis and SFD-BMD Generator

This project provides a MATLAB-based program and web application to perform structural analysis of simply supported Euler-Bernoulli beams subjected to point loads and uniformly distributed loads (UDLs). It calculates support reactions, shear force diagrams (SFD), bending moment diagrams (BMD), and beam deflection curves.

## Problem Statement

Analyze stress in simply supported beams with point loads and UDLs. Input beam length, load details, and positions to generate shear force, bending moment, and deflection diagrams for structural analysis and design validation.

## Objectives

- Create an interactive interface to input beam span, point loads, and UDL parameters.
- Compute support reactions using equilibrium equations.
- Calculate shear forces and bending moments along the beam.
- Obtain beam deflections using Euler-Bernoulli beam theory and numerical integration.
- Visualize SFD, BMD, and deflection curves with clear plots.
- Develop a user-friendly web application interface for easy access and visualization.

## Methodology

- Beam modeled as simply supported with loads including point forces and UDLs.
- Support reactions computed symbolically with equilibrium equations.
- Shear force and bending moment calculated by discretizing the beam and summing load effects.
- Deflection derived by numerical double integration of bending moment over EI (modulus of elasticity × moment of inertia) enforcing boundary conditions.
- MATLAB used for numerical calculations and plotting.
- Web application developed to allow online interaction with beam analysis.

## Usage

### MATLAB Program

- Input beam length, number of point loads, their magnitudes, and positions.
- Input number and details of UDLs (intensity and start-end positions).
- Run program to get textual results and visual diagrams for SFD, BMD, and deflection.
- Program outputs maximum deflection and locations of interest.

### Web Application

- Visit: https://beam-analysis-calculator.vercel.app
- Input beam parameters and loads interactively through the web interface.
- Visualize beam behavior and numerical results immediately.
- Suitable for quick design assessments and educational purposes.

## SFD-BMD Generator Script

- MATLAB script to compute shear force and bending moment diagrams for generic simply supported beams.
- Accepts beam span, multiple point loads, and multiple UDLs.
- Calculates reactions symbolically.
- Discretizes beam and integrates load effects.
- Produces plots marking key points such as zero shear, maximum bending moments, load locations and distributed load extents.

## Limitations

- Supports only simply supported beams with vertical point loads and uniform distributed loads.
- Does not currently handle other beam supports (cantilever, fixed).
- Deflection accuracy depends on numerical integration resolution.
- Future work may include GUI enhancements, other beam types, and non-uniform load types.

## Learning Outcomes

- Applied static equilibrium to determine support reactions.
- Used numerical methods for shear force and bending moment calculations.
- Integrated beam bending theory to compute deflections.
- Developed MATLAB skills for structural analysis.
- Created web-based tool for interactive beam visualization.

## References

- R. C. Hibbeler, Structural Analysis, 10th Edition, Pearson Education, 2016.
- S. P. Timoshenko, Strength of Materials, Part II, 3rd Edition, 1956.
- J. M. Gere and S. P. Timoshenko, Mechanics of Materials, 4th Edition, 1997.

---

Developed by Chandra Shekhar at Indian Institute of Technology Gandhinagar.

