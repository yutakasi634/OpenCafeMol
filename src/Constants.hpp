#ifndef OPEN_AICG2_PLUS_CONSTANTS_HPP
#define OPEN_AICG2_PLUS_CONSTANTS_HPP

#include "EnumTables.hpp"

// define physical constnt
static const double Na = 6.02214076e23;
static const double kB = 1.380649e-23; // J/K
static const double pi = 3.1415926535897932385;

// prepare fundamental constant for calculation
static const double cafetime = std::sqrt(1.0/OpenMM::KJPerKcal)*0.1; // Ps

// prepare parameter table for cubic spline of flexible local angle
static const std::map<AAType, std:array<double, 10>> =
    {
        {AAtype::ALA, {5.00,    1.34,    0.84,    1.17,    0.82,    1.00,    1.27,    1.52,    3.20,   10.00},
        {AAtype::ARG, {5.00,    1.43,    0.84,    1.09,    0.88,    1.01,    1.11,    1.60,    3.26,   10.00},
        {AAtype::ASN, {5.00,    1.38,    0.64,    0.92,    1.06,    1.27,    1.39,    1.71,    3.18,   10.00},
        {AAtype::ASP, {5.00,    1.32,    0.64,    0.94,    1.00,    1.24,    1.61,    1.81,    3.62,   10.00},
        {AAtype::CYS, {5.00,    1.59,    0.96,    1.00,    0.83,    0.97,    1.11,    1.44,    3.34,   10.00},
        {AAtype::GLN, {5.00,    1.41,    0.79,    1.05,    0.96,    1.01,    1.13,    1.66,    3.21,   10.00},
        {AAtype::GLU, {5.00,    1.30,    0.77,    1.09,    0.91,    1.05,    1.26,    1.76,    3.39,   10.00},
        {AAtype::GLY, {5.00,    1.60,    0.56,    1.21,    1.38,    1.12,    1.08,    1.23,    2.52,   10.00},
        {AAtype::HIS, {5.00,    1.50,    0.83,    0.99,    0.98,    1.05,    1.07,    1.43,    2.93,   10.00},
        {AAtype::ILE, {5.00,    1.79,    1.17,    0.92,    0.65,    0.89,    1.18,    2.06,    3.58,   10.00},
        {AAtype::LEU, {5.00,    1.54,    0.87,    1.07,    0.81,    0.88,    1.23,    2.11,    3.65,   10.00},
        {AAtype::LYS, {5.00,    1.36,    0.80,    1.04,    0.93,    1.02,    1.18,    1.72,    3.46,   10.00},
        {AAtype::MET, {5.00,    1.58,    0.87,    1.06,    0.89,    0.91,    1.09,    1.60,    4.12,   10.00},
        {AAtype::PHE, {5.00,    1.66,    0.92,    0.99,    0.91,    0.97,    1.01,    1.39,    3.23,   10.00},
        {AAtype::PRO, {5.00,    1.24,    0.82,    1.30,    0.64,    0.89,    2.23,    3.43,    4.79,   10.00},
        {AAtype::SER, {5.00,    1.39,    0.84,    1.21,    0.95,    0.99,    1.05,    1.29,    2.89,   10.00},
        {AAtype::THR, {5.00,    1.59,    0.84,    1.01,    0.99,    1.04,    0.96,    1.44,    3.29,   10.00},
        {AAtype::TRP, {5.00,    1.46,    0.84,    1.04,    0.88,    1.03,    1.14,    1.61,    3.25,   10.00},
        {AAtype::TYR, {5.00,    1.71,    0.91,    1.02,    0.90,    0.98,    0.98,    1.43,    3.26,   10.00},
        {AAtype::VAL, {5.00,    1.75,    1.19,    0.95,    0.72,    0.87,    1.01,    1.91,    3.59,   10.00},
    };

#endif // OPEN_AICG2PLUS_CONSTANTS_HPP
