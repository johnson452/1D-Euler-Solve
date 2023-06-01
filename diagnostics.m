function [] = diagnostics(rho,u,p,E,grid)


% Run only at select iterations:
if (mod ( grid.iter, grid.Output_interval ) == 0 || grid.iter == grid.NT)

    % Clear the figure
    clf()
    %Plot the diagnostic output comparison to fig1 GA Sod
    subplot(2,3,1)
    plot(grid.x,rho)
    title("Density")

    subplot(2,3,2)
    plot(grid.x,u)
    title("Velocity")

    subplot(2,3,3)
    plot(grid.x,p)
    title("Pressure")

    subplot(2,3,4)
    plot(grid.x,E - (1/2)*(u.*u))
    title("Energy")

    pause(0.01)

end



end