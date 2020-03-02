#!/usr/bin/python

# %%
import numpy
import meshio
mesh = meshio.read('../data/Ex11vec_noscale0385.vtk')
#print(vars(mesh))
points = mesh.points
#print(points)
pressure = mesh.point_data['pressure']
#pressure
velocity = mesh.point_data['Velocity']
#print(velocity[:,0])
#velocity

import plotly.io as pio
pio.renderers.default = "iframe"
import plotly.graph_objects as go
fig= go.Figure(data=go.Isosurface(
    x= points[:,0],
    z=points[:,1],
    y=points[:,2],
    colorscale=[[1,'rgb(0,255,0)'],[1,'rgb(255,0,0)']],
    showscale=False,
    lighting=dict(specular=0.5,ambient=0.8),
    value=velocity[:,0],
    isomin=-0.006,
    isomax=0.006,
    caps=dict(x_show=False, y_show=False,z_show=False)
))

fig.update_layout(
    title = dict(text="Isosurfaces of the most amplified streamwise velocity fluctuations", font = dict(family="CMU Serif", size=24,color="#000")),
    scene=dict(
      #  annotations=dict(
      #      x = 0,
      #      y = 0,
      #      z = 0,
       #     text="<i>y</i>"),
        camera=dict(
            eye=dict(x=1.25,y= 1,z = 1.25),
            projection=dict(type="orthographic")),
    xaxis = dict(
        tickangle=0,
        title=dict(
        text ="<i>x</i>",
        font=dict(
        family="CMU Serif",
        size = 20,color="#000"
    )),
    
    tickfont = dict(family="CMU Serif",size = 12,color="#000")),
    yaxis = dict(
        tickangle=0,
        title=dict(
        text ="<i>z</i>",
        font=dict(
        family="CMU Serif",
        size = 20,color="#000"
    )),
    tickfont = dict(family="CMU Serif",size = 12,color="#000")),
    zaxis = dict(
        tickangle=0,
        title=dict(
        text ="<i>y</i>",
        font=dict(
        family="CMU Serif",
        size = 20,color="#000"
    )),
    tickfont = dict(family='"CMU Serif", "Times New Roman",Arial ',size = 12,color="#000"))
))
fig.show()


# %%

# %%
