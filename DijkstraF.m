%% Group 02's Dijkstra's Algorithm
% This program will take an array of points which make up the map of ND and
% will tell the user the shortest possible distance from one point to
% another, and will also plot the results.
 
% Authors: Section 09, Group 02
% Date: 11 April 2017
 
% Referenced: dijkstra.m by Joseph Kirk
 
%% Set Function
 
function [dist,path] = DijkstraF(intersections, paths,startID,finishID)
%% Initialize parameters
% creates array of id's
intersectionIDs = intersections(:,1);
% find number of intersections and number of columns
[numIntersections,cols] = size(intersections);
% im not sure what sparse is if we're being real
table = sparse(num_map_pts,2);
% set all of the distances to infinity so shorter distances can be saved
shortest_dist = Inf(numIntersections,1);
% initizlize the array of points that have been checked
settled = zeros(numIntersections,1);
% Also not sure what this is
path = num2cell(NaN(numIntersections,1)); % oh initialize path output
% also not sure what this is
col = 2; % i would guess columns?
% also don't know what this is
pidx = find(startID == intersectionIDs);
% I assume this is initializing shortest distance
shortest_dist(pidx) = 0;
% initialze first column location maybe?
table(pidx,col) = 0;
% initialize first pidx location for settled array
settled(pidx) = 1;
% initilize pidx location for the path 
path(pidx) = {startID};
 
%% Calculate shortest path for all nodes
 
while sum(~settled) > 0
    % update table
    table(:,col-1) = table(:,col);
    table(pidx,col) = 0;
    % find id of nodes next to selected intersections
    neighborIDs = [paths(intersectionIDs(pidx) == paths(:,2),3);...
        paths(intersectionIDs(pidx) == paths(:,3),2)];
    % calculate distances to neighboring intersections
    for k = 1:length(intersectionIDs)
        % isolate the intersection is the same as the neighbor iteration
        cidx = find(neighborIDs(k) == intersectionIDs);
        if ~settled(cidx)
            % take the square root of the sum of the differences between
            % each of the intersection nodes that neighbor the selected
            % node squared to find the shortest one
            d = sqrt(sum((intersections(pidx,2:cols) -...
                intersections(cidx,2:cols)).^2));
            % checks if this intersection has been checked or if this point
            % that is being checked is a longer than the previous length +
            % the new distance calculated
            if (table(cidx,col-1) == 0) || ...
                    (table(cidx,col-1) > (table(pidx,col-1) + d))
                table(cidx,col) = table(pidx,col-1) + d;
                % sets temporary path to check as the path of pidx
                tempPath = path(pidx);
                % sets cidx path cell as first cell as the temporary
                % path and the id of the neighboring intersection to be
                % checked
                path(cidx) = {[tempPath{1} neighborIDs(k)]};
            else
                % if the intersection was checked or was found to be less,
                % move to next point in the table to check
                table(cidx,col) = table(cidx,col-1);
            end
        end
    end
    % find the minimum non-zero value in the table
    nidx = find(table(:,col));
    ndx = find(table(nidx,col) == min(table(nidx,col)));
    
    % if no possible route is found (the minimum value in the table is
    % still zero/empty, the code ends) (this shouldn't happen)
    if isempty(ndx)
        break
    else
        % sets the path list pidx as the nidx list at the minimum
        % value in the list
        pidx = nidx(ndx(1));
        % sets the shortest distance of the pidx list as the value from
        % the table to finish this off
        shortest_dist(pidx) = table(pidx,col);
        % settles the last point and calls it a day
        settled(pidx) = 1;
    end
end
% finds the value of the intersection ID that is the finishID
zz = find(finishID == intersectionID);
% sets the shortest distance for the output
dist = shortest_dist(zz);
% sets the path for the output
path = path(zz);
% sets the path to be the first path cell
path = path{1};








