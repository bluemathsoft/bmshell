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

import * as bluemath from 'bluemath'
import {plot} from './plot'

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

(<any>window).log = function (...args:any[]) {
  let s = '';
  for(let arg of args) {
    s += arg.toString().replace(/\n/g,'<br/>').replace(/\s/g,'&nbsp;');
  }
  let el = document.createElement('p');
  el.innerHTML = s;
  document.body.appendChild(el);
};


(<any>window).plot = plot;
