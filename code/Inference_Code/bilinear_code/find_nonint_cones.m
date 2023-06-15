function cone_set = find_nonint_cones(hyper_norms)

% Finds all non-intersecting cones created from the set of hyperplanes
% whos normals are input (all hyperplanes must include the origin).
% Algorithm starts with two hyperplanes, creating an initial 4 cones. Then
% each hyperplane is compared to each existing cone. If the hyperplane
% intersects the cone, then that cone is replaced by two new cones in the
% cone set. Otherwise, that cone is left alone. 
% 
% The algorithm uses Farkas Lemma at each step to determine if the
% intersection exists. Farkas Lemma states that only one of the following
% is true:
%
% 1) there exists an x s.t. A^Tx = b and x >= 0
% 2) there exists a y s.y. Ay >= 0 and b^ty < 0
%
% To show that there exists an [A; -b^T]y > 0 we show that there does not
% exist an x s.t. A^Tx = -b and x >= 0. 

% Normalize rows
hyper_norms = diag(1./sqrt(sum(hyper_norms.^2,2)))*hyper_norms;

% Set up initial cones (negatives are automatically included)
cone_set = {hyper_norms(1:2,:), [hyper_norms(1,:);-hyper_norms(2,:)]};


% Iteratively build up cone set by testing if and where each hyperplane splits existing cones
for kk = 3:size(hyper_norms,1)
    for ll = 1:numel(cone_set)
        % Test if hyperplane splits this cone
        n = size(cone_set{ll}, 1);
        cvx_begin quiet
	        variable x(n);
	        minimize( norm(cone_set{ll}.'*x - hyper_norms(kk,:).'));
	        subject to
	            x <= 0;
	    cvx_end
        if norm(cone_set{ll}.'*x - hyper_norms(kk,:).') >= 1e-8
	        % Condition is not violated, test the next
            cvx_begin quiet
	            variable x(n);
	            minimize( norm(cone_set{ll}.'*x + hyper_norms(kk,:).'));
	            subject to
	                x >= 0;
	        cvx_end
            if norm(cone_set{ll}.'*x + hyper_norms(kk,:).') >= 1e-8
            % Condition is still violated?
		        split_cone = 1;
            else
                split_cone = 0;
            end
        else
	        split_cone = 0;
        end
	if split_cone == 0
            % Don't split the cone
	else
            % Split the cone into two cones
            cone_set = {cone_set{1:ll-1},  [cone_set{ll};hyper_norms(kk,:)], [cone_set{ll};-hyper_norms(kk,:)], cone_set{ll+1:end}};
	end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%