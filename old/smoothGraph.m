% clc
% clear all
% skel = StdIP.readImg3D(1,1);
load skel.mat
% skeleton_graph = StdIP.skel2Graph3D(skel);
% 
% %%
% [S,C] = graphconncomp(skeleton_graph.SparseGraph,'Directed',false);
% all_node_lst = skeleton_graph.NodeList;
% all_parent_lst = skeleton_graph.ParentList;
% all_node_pos = skeleton_graph.NodePosition;
% %% Only for this component
% this_comp_id = find(C == 1);
% comp_node_lst = all_node_lst(this_comp_id);
% comp_parent_lst = skeleton_graph.ParentList(this_comp_id);
% comp_node_pos = skeleton_graph.NodePosition(this_comp_id);
% 
% %% Graph for this component
% N = length(this_comp_id);
% this_comp_graph = zeros(size(N));
% for ii = 1 : N
%     node_id = ii;
%     parent = comp_parent_lst(ii);
%     if parent == -1
%         parent_id(ii) = -1;
%     else
%         parent_id(ii) = find(comp_node_lst == parent);
%         this_comp_graph(ii,parent_id(ii)) = 1;
%         this_comp_graph(parent_id(ii),ii) = 1;
%     end
% end
% %% Lets spit out one swc for each component
% renamed_comp_node_id = 1:N;
% renamed_comp_parent_id = parent_id;
% renamed_comp_node_pos = comp_node_pos;
% %%
% [x,y,z] = ind2sub(size(skel),renamed_comp_node_pos);
% nodePos.x = x';
% nodePos.y = y';
% nodePos.z = z';
% 
% [splinedG,splinedNode] = graphSpline3D(sparse(this_comp_graph),nodePos,4);
%%
% splined_graph.SparseGraph = splinedG;
% N = size(splinedG,1);
% splined_graph.NodeList = 1:N;
% splined_graph.NodePosition = sub2ind(size(skel),splinedNode.x,splinedNode.y,splinedNode.z);

num_preserve = 3;
CC = bwconncomp(skel);
N = CC.NumObjects;
if N > num_preserve
    m = num_preserve;
else
    m = N;
end
for ii = 1 : m
    sub_skel = zeros(size(skel));
    skel_id = CC.PixelIdxList{ii};
    sub_skel(skel_id) = 1;
    skeleton_graph = StdIP.skel2Graph3D(sub_skel);
    [x,y,z] = ind2sub(size(skel),skeleton_graph.NodePosition);
    nodePos.x = y';
    nodePos.y = x';
    nodePos.z = z';
    
    [splinedG,splinedNode] = graphSpline3D(skeleton_graph.SparseGraph,nodePos,4);
    
    toSWC(skel,sparse(splinedG),splinedNode,6,2,2,1);
end


