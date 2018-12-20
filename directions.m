function [directions] = directions(pointvec,pointmat)
%DISTS computes the unit vector between one vector and a matrix of points
%   DISTS returns a matrix of unitary direction vectors corresponding to
%   the directions between coordinate vector POINTVEC and a matrix of
%   cooordinates POINTMAT
%
% Thomas Blum, Geo-Energy Lab, EPFL, October 2018

tmp = bsxfun(@minus,pointvec,pointmat);
directions = tmp./vecnorm(tmp);