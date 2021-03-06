#include <Spirit/Quantities.h>
#include <data/State.hpp>
#include <engine/Vectormath.hpp>

void Quantity_Get_Magnetization(State * state,  float m[3], int idx_image, int idx_chain)
{
    std::shared_ptr<Data::Spin_System> image;
    std::shared_ptr<Data::Spin_System_Chain> chain;
    from_indices(state, idx_image, idx_chain, image, chain);

	// image->Lock();
    auto mag = Engine::Vectormath::Magnetization(*image->spins);
	// image->Unlock();
	image->M = Vector3{ mag[0], mag[1], mag[2] };

    for (int i=0; i<3; ++i) m[i] = (float)mag[i];
}

float Quantity_Get_Topological_Charge(State * state, int idx_image, int idx_chain)
{
    std::shared_ptr<Data::Spin_System> image;
    std::shared_ptr<Data::Spin_System_Chain> chain;
    from_indices(state, idx_image, idx_chain, image, chain);

	// image->Lock();
    scalar charge = Engine::Vectormath::TopologicalCharge(*image->spins);
	// image->Unlock();

    return (float)charge;
}