    bool secondOrder = true;


// read Riemann solver coeffs
if(mesh.solutionDict().found("Riemann"))
{
    // Read Riemann coeffs dict
    dictionary riemann = mesh.solutionDict().subDict("Riemann");
    if (riemann.found("secondOrder"))
    {
        secondOrder = Switch(riemann.lookup("secondOrder"));
    }
}
