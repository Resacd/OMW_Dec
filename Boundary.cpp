#include "Boundary.h"
#include <cmath>
#include <iostream>

Boundary::Boundary(int coord_x, int coord_y, int coord_z, std::shared_ptr <VelocitySet> velset, NodeType nt)
    : Node(coord_x, coord_y, coord_z, nt) 
{

    m_inletVel.resize(_dim);
    velSet = velset;

}

Boundary::~Boundary() {


}

double Boundary::Density() {

	return 1;
}


void Boundary::SetInletVel(double ux, double uy, double uz)
{
    m_inletVel.push_back(ux);
    m_inletVel.push_back(uy);
    m_inletVel.push_back(uz);
}

std::vector<double> Boundary::Velocity(std::shared_ptr<VelocitySet> velSet)
{
    if (BoundaryType_ == BoundaryType::Inlet) {
        return m_inletVel;
    }
    if (BoundaryType_ == BoundaryType::Wall) {

        return { 0.0,0.0,0.0 };
    }
}

void Boundary::ApplyEqBc()
{
    for (int dir = 0; dir < _nLatNodes; dir++) {
        double rho = 1.0;
        double du = m_inletVel[0] * velSet->GetDirection(dir)[0] +
            m_inletVel[1] * velSet->GetDirection(dir)[1] +
            m_inletVel[2] * velSet->GetDirection(dir)[2];
        double u_sqr = m_inletVel[0] * m_inletVel[0] + m_inletVel[1] * m_inletVel[1] + m_inletVel[2] * m_inletVel[2];
        double equilibriumDist = velSet->GetWeight(dir) * rho *
            (1 + 3 * du + 9.0 / 2.0 * du * du - 3.0 / 2.0 * u_sqr);

        m_distributions[dir] = equilibriumDist;
    }
}


void Boundary::MacroscopicValue(std::shared_ptr<VelocitySet> velSet) {



    density = 0;
    velocity[0] = 0.0;
    velocity[1] = 0.0;
    velocity[2] = 0.0;

    for (int dir = 0; dir < _nLatNodes; dir++) {

        density += m_distributions[dir];
        velocity[0] += m_distributions[dir] * velSet->GetDirection(dir)[0];
        velocity[1] += m_distributions[dir] * velSet->GetDirection(dir)[1];
        velocity[2] += m_distributions[dir] * velSet->GetDirection(dir)[2];

    }

    velocity[0] /= density;
    velocity[1] /= density;
    velocity[2] /= density;
    //density = 1;
    //velocity = m_inletVel;

}

#ifdef D3Q27

void Boundary::ApplyBc()
{
    double rho = ((m_distributions[0] + m_distributions[3] + m_distributions[4]       +
                   m_distributions[5] + m_distributions[6] + m_distributions[15]      +
                   m_distributions[16] + m_distributions[17] + m_distributions[18])   +
                   2 * (m_distributions[2] + m_distributions[9] + m_distributions[10] +
                   m_distributions[12] + m_distributions[14] + m_distributions[20]    +
                   m_distributions[22] + m_distributions[24] + m_distributions[26])) / (1 - m_inletVel[0]);

    double Ny = (1.0 / 6) * (m_distributions[3] + m_distributions[15] + m_distributions[16]
        - m_distributions[4] - m_distributions[17] - m_distributions[18]) - (1.0 / 9) * (rho * m_inletVel[1]);

    double Nz = (1.0 / 6) * (m_distributions[5] + m_distributions[15] + m_distributions[17]
        - m_distributions[6] - m_distributions[16] - m_distributions[18]) - (1.0 / 9) * (rho * m_inletVel[2]);


    m_distributions[1] = m_distributions[2]   + (4 / 9.0) * rho  * m_inletVel[0];

    m_distributions[7] = m_distributions[10]  + (1 / 9.0) * rho  * (m_inletVel[0] + m_inletVel[1]) - Ny;
                                                                    
    m_distributions[8] = m_distributions[9]   + (1 / 9.0) * rho  * (m_inletVel[0] - m_inletVel[1]) + Ny;
                                                                    
    m_distributions[11] = m_distributions[14] + (1 / 9.0) * rho  * (m_inletVel[0] + m_inletVel[2]) - Nz;
                                                                    
    m_distributions[13] = m_distributions[12] + (1 / 9.0) * rho  * (m_inletVel[0] - m_inletVel[2]) + Nz;
                                                                    
    m_distributions[19] = m_distributions[20] + (1 / 36.0) * rho * (m_inletVel[0] + m_inletVel[1] + m_inletVel[2]) - Ny - Nz;
                                                         
    m_distributions[21] = m_distributions[22] + (1 / 36.0) * rho * (m_inletVel[0] + m_inletVel[1] - m_inletVel[2]) - Ny + Nz;
                                                         
    m_distributions[23] = m_distributions[24] + (1 / 36.0) * rho * (m_inletVel[0] - m_inletVel[1] + m_inletVel[2]) + Ny - Nz;
                                                         
    m_distributions[25] = m_distributions[26] + (1 / 36.0) * rho * (m_inletVel[0] - m_inletVel[1] - m_inletVel[2]) + Ny + Nz;


    MacroscopicValue(velSet);
}

#endif

#ifdef D3Q19

void Boundary::ApplyBc(double& vel)
{

    m_inletVel[0] = 0.0;
    m_inletVel[1] = 0.0;
    m_inletVel[2] = vel;

    double rhoN = (m_distributions[0] + m_distributions[1] + m_distributions[4] + m_distributions[5] +
        m_distributions[8] + m_distributions[9] + m_distributions[12] + m_distributions[13] + m_distributions[18] +
        2 * (m_distributions[3] + m_distributions[7] + m_distributions[11] + m_distributions[16] + m_distributions[17]))/(1 - m_inletVel[1]);

    double Nx = (1.0 / 2) * (m_distributions[0] + m_distributions[8] + m_distributions[9]
        - m_distributions[1] - m_distributions[12] - m_distributions[13]) - (1.0 / 3) * (rhoN * m_inletVel[0]);

    double Nz = (1.0 / 2) * (m_distributions[4] + m_distributions[8] + m_distributions[12]
        - m_distributions[5] - m_distributions[9] - m_distributions[13]) - (1.0 / 3) * (rhoN * m_inletVel[2]);


  m_distributions[2] = m_distributions[3]   + (1 / 3.0) * rhoN * m_inletVel[1];

  m_distributions[6] = m_distributions[11]  + (1 / 6.0) * rhoN * (m_inletVel[0] + m_inletVel[1]) - Nx;

  m_distributions[10] = m_distributions[7]  + (1 / 6.0) * rhoN * (m_inletVel[1] - m_inletVel[0]) + Nx;

  m_distributions[14] = m_distributions[17]  + (1 / 6.0) * rhoN * (m_inletVel[1] + m_inletVel[2]) - Nz;

  m_distributions[15] = m_distributions[16]  + (1 / 6.0) * rhoN * (m_inletVel[1] - m_inletVel[2]) + Nz;

 // MacroscopicValue(velSet);
}

#endif


