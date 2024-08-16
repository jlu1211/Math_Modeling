% Define the system of differential equations
f = @(t, Y, c) [Y(1) * (3 - Y(1) - 2 * Y(2)) + c * Y(1);
                Y(2) * (2 - Y(1) - Y(2)) + c * Y(2)];

% Fixed points labels
labels = {'Node (stable)', 'Node (unstable)', 'Saddle', 'Spiral (unstable)'};

% Values of c to be used
c_values = [-10, -4-2*sqrt(2), -5, -3, -2.5, -2, -1.5, -4+2*sqrt(2), -1.15, -1, 0];

% Plot the vector field for each fixed point with x and y range +- 0.5 around the fixed point
for c = c_values
    % Fixed points
    fixed_points = [0, 0; 0, 2 + c; 3 + c, 0; 1 + c, 1];
    
    for j = 1:size(fixed_points, 1)
        % Create a grid of (x, y) points around the fixed point
        x_range = fixed_points(j, 1) + (-0.5:0.05:0.5);
        y_range = fixed_points(j, 2) + (-0.5:0.05:0.5);
        [x, y] = meshgrid(x_range, y_range);

        % Initialize the vector field
        u = zeros(size(x));
        v = zeros(size(y));

        % Calculate the vectors at each point
        for i = 1:numel(x)
            Yprime = f(0, [x(i); y(i)], c);
            u(i) = Yprime(1);
            v(i) = Yprime(2);
        end

        % Plot the vector field with bigger arrows
        figure
        quiver(x, y, u, v, 5, 'b') % Increase the scaling factor to 1.5
        hold on
        plot(fixed_points(j, 1), fixed_points(j, 2), 'ro')

        % Annotate the fixed point
        text(fixed_points(j, 1), fixed_points(j, 2), labels{j}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')

        % Set plot limits and labels
        xlim([fixed_points(j, 1) - 0.5, fixed_points(j, 1) + 0.5])
        ylim([fixed_points(j, 2) - 0.5, fixed_points(j, 2) + 0.5])
        xlabel('x (Sheep Population)')
        ylabel('y (Rabbit Population)')
        title(['Phase Plane Analysis at Fixed Point (' num2str(fixed_points(j, 1)) ', ' num2str(fixed_points(j, 2)) '), c = ' num2str(c) ''])
        grid on
        hold off
    end
end
