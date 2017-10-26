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

import * as THREE from 'three'
import {OrbitControls} from 'three-orbitcontrols-ts'
import {TypedArray} from '@bluemath/common'

declare let Plotly:any; // Workaround to get rid of TypeScript compile error

export interface TessFormat3D {
  points? : Array<number[]>;
  line? : Array<number[]>;
  mesh? : {
    vertices : TypedArray|number[];
    faces : TypedArray|number[];
  },
  color? : number,
  targetGL? : number
}


function makeAxes() {
  var L = 50;
  var xAxisGeometry = new THREE.Geometry();
  var yAxisGeometry = new THREE.Geometry();
  var zAxisGeometry = new THREE.Geometry();
  xAxisGeometry.vertices.push(
    new THREE.Vector3(0,0,0),
    new THREE.Vector3(L,0,0),
    new THREE.Vector3(0,0,0),
    new THREE.Vector3(-L,0,0)
  );
  yAxisGeometry.vertices.push(
    new THREE.Vector3(0,0,0),
    new THREE.Vector3(0,L,0)
  );
  zAxisGeometry.vertices.push(
    new THREE.Vector3(0,0,0),
    new THREE.Vector3(0,0,L),
    new THREE.Vector3(0,0,0),
    new THREE.Vector3(0,0,-L)
  );
  xAxisGeometry.computeLineDistances();
  yAxisGeometry.computeLineDistances();
  zAxisGeometry.computeLineDistances();

  return [
    new THREE.LineSegments(xAxisGeometry,
      new THREE.LineDashedMaterial({
        color: 0xffa805, dashSize: 0.3, gapSize: 0.15, linewidth: 1
      })),
    new THREE.LineSegments(yAxisGeometry,
      new THREE.LineDashedMaterial({
        color: 0xff05b4, dashSize: 0.3, gapSize: 0.15, linewidth: 1
      })),
    new THREE.LineSegments(zAxisGeometry,
      new THREE.LineDashedMaterial({
        color: 0x05ff9b, dashSize: 0.3, gapSize: 0.15, linewidth: 1
      }))
  ]
}

interface GL {
  scene : THREE.Scene;
  canvas : HTMLCanvasElement;
  glrndr : THREE.WebGLRenderer;
  camera : THREE.Camera;
  orbitControls : OrbitControls;
  cpointSprite : THREE.Texture;
}

function makeGL(canvas3:HTMLCanvasElement,width:number,height:number) : GL {
  let renderer = new THREE.WebGLRenderer({canvas:canvas3, antialias:true});
  renderer.setSize(width, height);
  renderer.setPixelRatio(window.devicePixelRatio);
  renderer.setClearColor(0xffffff, 1);

  let camera = new THREE.PerspectiveCamera(75, width/height, 1, 200);

  let scene = new THREE.Scene();
  //scene.fog = new THREE.Fog(0xffffff, 150,200);

  makeAxes().forEach(function (o) { scene.add(o); });

  let ambientLight = new THREE.AmbientLight(0x111111);
  scene.add(ambientLight);

  let ypluslight = new THREE.PointLight( 0xFFFFFF );
  ypluslight.position.set(0, 10, 0);
  scene.add(ypluslight);
  let yminuslight = new THREE.PointLight( 0xFFFFFF );
  yminuslight.position.set(0, -10, 0);
  scene.add(yminuslight);
  let amblight = new THREE.AmbientLight(0x444444);
  scene.add(amblight);

  camera.position.x = 4;
  camera.position.y = 4;
  camera.position.z = 4;

  var orbitControls = new OrbitControls(camera, renderer.domElement);
  orbitControls.enableDamping = true;
  orbitControls.dampingFactor = 0.25;
  orbitControls.enableZoom = true;
  orbitControls.autoRotate = false;
  orbitControls.autoRotateSpeed = 0.35;

  return {
    cpointSprite : THREE.ImageUtils.loadTexture( "cpoint.png" ),
    scene : scene,
    canvas : canvas3,
    glrndr : renderer,
    camera : camera,
    orbitControls : orbitControls
  };
}

export class Renderer {

  private div : HTMLElement;
  private plotlyLayout : any;
  private width : number;
  private height : number;
  private gl1 : GL;
  private gl2? : GL;

  constructor(
    div:HTMLElement,
    type:'plotly'|'threejs',
    size={width:400,height:300},
    split=false)
  {
    this.div = div;
    this.width = size.width;
    this.height = size.height;
    div.style.width = size.width+'';
    div.style.height = size.height+'';

    if(type === 'plotly') {
      this.initPlot(split);
    } else {
      this.initGL(split);
    }
  }

  private initPlot(split:boolean) {
    let margin = 0.05;
    if(split) {
      this.plotlyLayout = {};
      if(this.width < this.height) { // Portrait - Vertical stack
        this.plotlyLayout.xaxis = { anchor:'y1' };
        this.plotlyLayout.xaxis2 = { anchor:'y2' };
        this.plotlyLayout.yaxis = { domain:[0.5+margin,1] };
        this.plotlyLayout.yaxis2 = { domain:[0,0.5-margin] };
      } else { // Landscape - Horizontal stack
        console.assert(false); // TODO
        this.plotlyLayout.xaxis = { domain:[0.5+margin,1] };
        this.plotlyLayout.xaxis2 = { domain:[0,0.5-margin] };
        this.plotlyLayout.yaxis = { };
        this.plotlyLayout.yaxis2 = { anchor : 'x2' };
      }
      this.plotlyLayout.margin = { t:0, b:0, l:0, r:0 };
      this.plotlyLayout.width = this.width;
      this.plotlyLayout.height = this.height;
    } else {
      this.plotlyLayout = {};
      this.plotlyLayout.xaxis = {anchor:'y1'};
      this.plotlyLayout.yaxis = {};
      this.plotlyLayout.margin = { t:20, b:20, l:20, r:20 };
      this.plotlyLayout.width = this.width;
      this.plotlyLayout.height = this.height;
      this.plotlyLayout.showlegend = false;
    }
  }

  private initGL(split:boolean) {
    if(split) {
      let canvas1 = document.createElement('canvas');
      this.div.appendChild(canvas1);
      this.gl1 = makeGL(canvas1,this.width,this.height/2);
      let canvas2 = document.createElement('canvas');
      this.div.appendChild(canvas2);
      this.gl2 = makeGL(canvas2,this.width,this.height/2);
    } else {
      let canvas1 = document.createElement('canvas');
      this.div.appendChild(canvas1);
      this.gl1 = makeGL(canvas1,this.width,this.height);
    }
  }

  render2D(traces:any,options?:any) {
    if(options && options.range) {
      this.plotlyLayout.xaxis.range = options.range.x;
      this.plotlyLayout.yaxis.range = options.range.y;
      if(this.plotlyLayout.xaxis2) {
        this.plotlyLayout.xaxis2.range = options.range.x;
      }
      if(this.plotlyLayout.yaxis2) {
        this.plotlyLayout.yaxis2.range = options.range.y;
      }
    }
    if(options && options.ismatrix) {
      this.plotlyLayout.yaxis.autorange = 'reversed';
      this.plotlyLayout.xaxis.side = 'top';
      this.plotlyLayout.xaxis.ticklen = 0;
      this.plotlyLayout.yaxis.ticklen = 0;
    }
    if(options && options.title) {
      this.plotlyLayout.title = options.title;
      this.plotlyLayout.margin.t = 40;
    }
    Plotly.newPlot(this.div, traces, this.plotlyLayout, {staticPlot:true});
  }

  render3D(arrtess:TessFormat3D[], options?:any) {
    for(let tess of arrtess) {
      this.render3DSingle(tess,options);
    }
  }

  render3DSingle(tess:TessFormat3D,options?:any) {

    let gl:GL;
    if(tess.targetGL && tess.targetGL === 2) {
      console.assert(!!this.gl2);
      gl = this.gl2!;
    } else {
      gl = this.gl1;
    }

    if(options && options.title) {
      let titleText = document.createElement('div');
      titleText.style.width = this.width+'px';
      titleText.style.textAlign = 'center';
      titleText.style.fontFamily = 'arial';
      titleText.style.marginBottom = '3px';
      titleText.style.marginTop = '3px';
      titleText.textContent = options.title;
      this.div.insertBefore(titleText,gl.canvas);
    }

    // Add geometry mesh
    if(tess.mesh) {
      let loader = new THREE.JSONLoader();
      let geometry = loader.parse(tess.mesh).geometry;
      var material = new THREE.MeshLambertMaterial({
        color: tess.color || 0x4444ff,
        side: THREE.DoubleSide,
        flatShading : false,
        transparent : true,
        opacity : 0.5
      });
      geometry.computeVertexNormals();
      let mesh = new THREE.Mesh(geometry, material);
      gl.scene.add(mesh);
    }

    // Add points as sprites
    if(tess.points) {
      let cpointGeometry = new THREE.Geometry();
      tess.points.forEach(function (p) {
        cpointGeometry.vertices.push(new THREE.Vector3(p[0],p[1],p[2]));
      });
      let cpointMaterial = new THREE.PointsMaterial({
        size: 10, sizeAttenuation: false,
        map: gl.cpointSprite, alphaTest: 0.5, transparent: true });
      cpointMaterial.color.setHSL(1.0, 0.3, 0.7);

      let cpointParticles = new THREE.Points(cpointGeometry, cpointMaterial);
      gl.scene.add(cpointParticles);
    }

    if(tess.line) {
      let lineGeom = new THREE.Geometry();
      for(let i=0; i<tess.line.length-1; i++) {
        let p = tess.line[i];
        let n = tess.line[i+1];
        lineGeom.vertices.push(new THREE.Vector3(p[0],p[1],p[2]));
        lineGeom.vertices.push(new THREE.Vector3(n[0],n[1],n[2]));
      }
      lineGeom.computeLineDistances();
      let lineseg = new THREE.LineSegments(lineGeom, new THREE.LineBasicMaterial({
        color : 0xff0000, linewidth:2
      }));
      gl.scene.add(lineseg);
    }

    let orbitControls = gl.orbitControls;
    let glrndr = gl.glrndr;
    let scene = gl.scene;
    let camera = gl.camera;

    var render = function () {
      orbitControls.update();
      glrndr.render(scene, camera);
      requestAnimationFrame(render);
    };

    render();
  }
}
