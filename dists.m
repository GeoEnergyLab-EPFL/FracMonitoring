function [dists] = dists(pointvec,pointmat)
%DISTS computes the distance between one vector and a matrix of points
%   DISTS returns a vector corresponding to the distances between one
%   coordinate vector POINTVEC and a matrix of cooordinates POINTMAT
%
% Thomas Blum, Geo-Energy Lab, EPFL, October 2018

tmp = bsxfun(@minus,pointvec,pointmat);
dists = cellfun(@norm,num2cell(tmp,1));