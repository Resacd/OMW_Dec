#include "Node.h"

Node::Node(int coord_x, int coord_y, int coord_z, NodeType nt)
	: x_position(coord_x), y_position(coord_y), z_position(coord_z), NodeType_(nt)
{
	NodeType_ = nt;
	m_eqdistributions.resize(_nLatNodes);
        velocity.resize(_dim);
}


Node::~Node()
{
}

void Node::Stream(int dir)
{

   m_neighbours[dir]->m_distributions[dir] = m_newDistributions[dir];
}

double Node::Equilibrium(std::shared_ptr<VelocitySet> velSet, int dir)
{
	double du = Velocity(velSet)[0] * velSet->GetDirection(dir)[0] + Velocity(velSet)[1] * velSet->GetDirection(dir)[1]+ Velocity(velSet)[2] * velSet->GetDirection(dir)[2];
	double u_sqr = Velocity(velSet)[0] * Velocity(velSet)[0] + Velocity(velSet)[1] * Velocity(velSet)[1] + Velocity(velSet)[2] * Velocity(velSet)[2];

	return velSet->GetWeight(dir) * Density() * (1 + 3 * du + 9 / 2 * du * du - 3 / 2 * u_sqr);
}

void Node::Initialize(std::shared_ptr<VelocitySet> velSet, int dir)
{
	double rho = 1.0;   

	m_distributions.push_back(velSet->GetWeight(dir));
	m_newDistributions.push_back(0.0);
	velocity.resize(3);
}

double Node::Density()
{
	double density = 0;
	for (int dir = 0; dir < _nLatNodes; dir++)
	{
		density += this->m_distributions[dir];	
	}
	return density;
}

//std::vector<double> Node::Velocity(std::shared_ptr<VelocitySet> velSet)
//{
//	std::vector<double> velocity(_dim, 0.0);
//
//	for (int dim = 0; dim < _dim; dim++)
//	{
//
//		for (int dir = 0; dir < _nLatNodes; dir++)
//		{
//			velocity[dim] += velSet->GetDirection(dir)[dim] * this->m_distributions[dir];
//		}
//
//		velocity[dim] /= Density();
//	}
//
//	return velocity;
//}


std::vector<double> Node::Velocity(std::shared_ptr<VelocitySet> velSet)
{

	for (int dir = 0; dir < _nLatNodes; dir++)
	{
		for (int dim = 0; dim < 3; dim++)
		{
			velocity[dim] += velSet->GetDirection(dir)[dim] * this->m_distributions[dir];
		}
	}

	for (int dim = 0; dim < 3; dim++)
	{
		velocity[dim] /= Density();
	}

	return velocity;
}



void Node::MacroscopicValue(std::shared_ptr<VelocitySet> velSet) {

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
}



double Node::Pressure(std::shared_ptr<VelocitySet> velSet)
{
	double density = Density();
	double c_square = velSet->Get_c_square(0);
	double pressure = c_square * density;
	return pressure;
}


void Node::BounceBack(int dir, int opp_dir)
{

	m_distributions[opp_dir] = m_newDistributions[dir];
}

void Node::ApplyEqBc()
{

}

void Node::ApplyBc(double &vel)
{

}
