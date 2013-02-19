clc; clear all; close all;
%% initialize
%1 well of a well plate
nc = 10000; %num cells

%TIME PARAMS
dt = 1; %seconds
tend_h = 24;

tend_m = tend_h*60; tend_s = tend_m*60;
tv = 0:dt:tend_s; %time vector

%UBIQUITOUS CELL INFO
%RADII: assume spherical, equi-radii
r = 20.8153e-6; %cell radius, m
d = 2*r;
rad = num2cell(r*ones(1,nc)); %all cells have equal radius

%SEED CELLS
c
% %POSITION DATA: rectangular grid

%TNF alpha concentration