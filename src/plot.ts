/*

Copyright (C) 2017 Jayesh Salvi, Blue Math Software Inc.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*/

import * as bluemath from 'bluemath';
import {Renderer} from './renderer'
import {GeometryAdapter} from './adapter'
let {NDArray} = bluemath.common;

export interface PlotSpec {
  data? : any;
  type? : 'line'|'points'|'matrix';
  width? : number;
  height? : number;
  title? : string;
  interactive? : boolean;
}

export type PlotData = Array<any> |
  bluemath.common.NDArray |
  bluemath.geom.nurbs.BezierCurve |
  bluemath.geom.nurbs.BSplineCurve |
  bluemath.geom.nurbs.BezierSurface |
  bluemath.geom.nurbs.BSplineSurface;

function plot2DArray(object:bluemath.common.NDArray, spec?:PlotSpec) {
  let type:string = 'scatter';
  let mode:string|undefined = 'lines';
  if(spec && spec.type) {
    if(spec.type === 'points') {
      type = 'scatter';
      mode = 'markers';
    } else if(spec.type === 'line') {
      type = 'scatter';
      mode = 'lines';
    } else if(spec.type === 'matrix') {
      type = 'heatmap';
      mode = undefined;
    }
  }

  if(type === 'scatter' && object.shape[1] !== 2) {
    throw new Error('Data is not 2D point set');
  }

  let width = 300;
  let height = 200;
  if(spec && spec.width) {
    width = spec.width;
  }
  if(spec && spec.height) {
    height = spec.height;
  }

  let div = document.createElement('div');
  div.classList.add('bm-plot');
  document.body.appendChild(div);
  let rndr = new Renderer(div,'plotly',{width,height});
  let trace;
  if(type === 'scatter') {
    trace = {
      x : Array.from(object.getA(':',0).data),
      y : Array.from(object.getA(':',1).data),
      xaxis : 'x1',
      yaxis : 'y1',
      type : type,
      mode : mode
    };
    rndr.render2D([trace], {
      title : spec ? spec.title : undefined
    });
  } else if(type === 'heatmap') {
    trace = {
      z : object.toArray(),
      type : type
    };
    rndr.render2D([trace],{
      ismatrix : true,
      title : spec ? spec.title : undefined
    });
  }
  return div;
}

type GeometryType = bluemath.geom.nurbs.BezierCurve |
    bluemath.geom.nurbs.BSplineCurve |
    bluemath.geom.nurbs.BezierSurface |
    bluemath.geom.nurbs.BSplineSurface;

function plotGeometry(geometry:GeometryType, spec?:PlotSpec) {
  let div = document.createElement('div');
  div.classList.add('bm-plot');
  let width = (spec && spec.width) || 300;
  let height = (spec && spec.height) || 200;
  let rndr = new Renderer(
    div, is3D(geometry) ? 'threejs' : 'plotly',{width,height});
  new GeometryAdapter(geometry,rndr,{
    title : spec ? spec.title : undefined
  });
  return div;
}

function isGeometry(object:any) {
  return (
    object instanceof bluemath.geom.nurbs.BezierCurve ||
    object instanceof bluemath.geom.nurbs.BSplineCurve ||
    object instanceof bluemath.geom.nurbs.BezierSurface ||
    object instanceof bluemath.geom.nurbs.BSplineSurface
  );
}

function is3D(object:PlotData) {
  if(isGeometry(object)) {
    return (<any>object).dimension === 3;
  } else if(Array.isArray(object)) {
    return Array.isArray(object[0]) && object[0].length === 3;
  } else if(object instanceof NDArray) {
    return (<bluemath.common.NDArray>object).shape.length === 2 &&
      (<bluemath.common.NDArray>object).shape[1] === 3;
  }
  return false;
}

export function plot(data:PlotData|PlotData[], spec?:PlotSpec|PlotSpec[]) {
  if(!Array.isArray(data)) {
    data = [data];
  }
  if(spec && !Array.isArray(spec)) {
    spec = [spec];
  }

  let singleSpec = null;
  if(spec) {
    if(data.length === 1 && spec.length === 1) {
      singleSpec = spec[0];
    } else if(data.length === spec.length) {
      singleSpec = null;
    } else {
      throw new Error('Mismatch in data and spec count');
    }
  }

  let groupDiv = document.createElement('div');
  groupDiv.classList.add('bm-plot-group');
  document.body.appendChild(groupDiv);

  for(let i=0; i<data.length; i++) {
    let div;
    let plotspec = singleSpec ? singleSpec : (spec ? spec[i] : undefined);
    if(data[i] instanceof NDArray) {
      div = plot2DArray(data[i], plotspec);
    } else if(isGeometry(data[i])) {
      div = plotGeometry(data[i], plotspec);
    } else {
      console.warn('Plotting not supported for '+data[i]);
      continue;
    }
    groupDiv.appendChild(div);
  }

}
