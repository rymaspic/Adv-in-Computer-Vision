%% demo

clc;
clear;

source = imread('texture.jpg');       

target = SynthTexture(source, 13, [100, 100]);

imshow(target);