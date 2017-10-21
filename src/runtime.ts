
import * as bluemath from 'bluemath'
import {Renderer} from './renderer'
import {GeometryAdapter} from './adapter'

(<any>window).bluemath = bluemath;

(<any>window).console = {
  log : function(...args:any[]) {
    let msg:string = args.map(a => new String(a)).join('');
    let el = document.createElement('p');
    el.textContent = msg;
    document.body.appendChild(el);
  },
  error : function(...args:any[]) {
    let msg:string = args.map(a => new String(a)).join('');
    let el = document.createElement('p');
    el.style.color = '#ff0000';
    el.textContent = msg;
    document.body.appendChild(el);
  },
  warn : function(...args:any[]) {
    let msg:string = args.map(a => new String(a)).join('');
    let el = document.createElement('p');
    el.style.color = '#ff7700';
    el.textContent = msg;
    document.body.appendChild(el);
  },
  assert : function(condition:boolean) {
    if(!condition) {
      throw new Error("Assertion failed");
    }
  },
  clear : function () {
    let children = document.body.children;
    for(let i=children.length-1; i>=0; i--) { children[i].remove(); }
  }
};

(<any>window).scopedEval = function (codestr:string) {
  try {
    console.clear();
    eval(codestr);
  } catch(e) {
    console.error(e.toString());
  }
};

(<any>window).bmlog = function (...args:any[]) {
  let s = '';
  for(let arg of args) {
    s += arg.toString().replace(/\n/g,'<br/>').replace(/\s/g,'&nbsp;');
  }
  let el = document.createElement('p');
  el.innerHTML = s;
  document.body.appendChild(el);
};

function plot2DArray(object:bluemath.common.NDArray, options?:any) {
  let type:string = 'scatter';
  let mode:string|undefined = 'lines';
  if(options && options.type) {
    if(options.type === 'points') {
      type = 'scatter';
      mode = 'markers';
    } else if(options.type === 'line') {
      type = 'scatter';
      mode = 'lines';
    } else if(options.type === 'matrix') {
      type = 'heatmap';
      mode = undefined;
    }
  }

  if(type === 'scatter' && object.shape[1] !== 2) {
    throw new Error('Data is not 2D point set');
  }

  let width = 300;
  let height = 200;
  if(options && options.width) {
    width = options.width;
  }
  if(options && options.height) {
    height = options.height;
  }

  let div = document.createElement('div');
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
    rndr.render2D([trace]);
  } else if(type === 'heatmap') {
    trace = {
      z : object.toArray(),
      type : type
    };
    rndr.render2D([trace],{reverseYRange:true});
  }
}

(<any>window).bmplot = function (...args:any[]) {
  for(let i=0; i<args.length; i+=2) {
    let object = args[i];
    let options:any = (i < args.length-1) ? args[i+1] : {};
    if(object instanceof bluemath.common.NDArray) {
      let arr = <bluemath.common.NDArray>object;
      if(arr.shape.length === 2) {
        plot2DArray(object,options);
      } else {
        console.warn('Unplottable object '+object);
      }
    } else if(
      object instanceof bluemath.geom.nurbs.BezierCurve ||
      object instanceof bluemath.geom.nurbs.BSplineCurve ||
      object instanceof bluemath.geom.nurbs.BezierSurface ||
      object instanceof bluemath.geom.nurbs.BSplineSurface)
    {
      let div = document.createElement('div');
      document.body.appendChild(div);
      let width = (options && options.width) || 300;
      let height = (options && options.height) || 200;
      let rndr = new Renderer(div,'plotly',{width,height});
      new GeometryAdapter(object,rndr);
    } else {
      console.warn('Unplottable object '+object);
    }
  }
};
