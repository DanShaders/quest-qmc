#include "SRC/sema/FreeformGeometry.h"

namespace dqmc::sema {

namespace {

struct BipartitionBuildingContext : Context {
    BipartitionBuildingContext(Context& ctx)
        : Context(ctx)
    {
    }

    bool dfs(int u, int color);
    void find_bipartition();

    int n = lattice.sites.size();
    std::vector<int> partition = std::vector(n, 0);
};

// Colors a connected component of the lattice.
bool BipartitionBuildingContext::dfs(int u, int color)
{
    VERIFY(partition[u] == 0);

    partition[u] = color;

    for (int v = 0; v < n; ++v) {
        if (hamiltonian.hoppings[0](u, v) != 0 || hamiltonian.hoppings[1](u, v) != 0) {
            if (partition[v] == color) {
                return false;
            } else if (partition[v] == 0 && !dfs(v, -color)) {
                return false;
            }
        }
    }

    return true;
}

// Tries to color all connected components of the lattice, and if succeeds, saves the bipartition.
void BipartitionBuildingContext::find_bipartition()
{
    for (int i = 0; i < n; ++i) {
        if (partition[i] == 0 && !dfs(i, 1)) {
            return;
        }
    }

    bipartition.emplace(std::move(partition));
}

} // namespace

void find_bipartition_if_exists(Context& ctx)
{
    BipartitionBuildingContext { ctx }.find_bipartition();
}

} // namespace dqmc::sema
