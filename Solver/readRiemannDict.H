if(mesh.solutionDict().found("Riemann"))
{
    dictionary riemann = mesh.solutionDict().subDict("Riemann");
    if (riemann.found("secondOrder"))
    {
        secondOrder = Switch(riemann.lookup("secondOrder"));
    }
}
