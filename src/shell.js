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

let W = window.innerWidth-10;
let H = window.innerHeight-10;
var editor;
let containerElem = document.getElementById('container');
let editorElem = document.getElementById('editor');
let runtimeElem = document.getElementById('runtime');
let docholderElem = document.getElementById('docholder');
let docsCloseElem = document.getElementById('docsclose');
let exmChoices = document.getElementById('examples-choices');
let clearElem = document.getElementById('btn-clear');
let saveElem = document.getElementById('btn-save');

function nameToKey(name) {
  return name.replace(/[\(\)\s]+/g,'-').toLowerCase();
}

let preambleTypeDeclaration =
`
declare function log(...args:any[]) {}
declare interface PlotSpec {
  data? : any;
  type? : 'line'|'points'|'matrix';
  width? : number;
  height? : number;
  title? : string;
  interactive? : boolean;
}

declare type PlotData = Array<any> |
  bluemath.common.NDArray |
  bluemath.geom.nurbs.BezierCurve |
  bluemath.geom.nurbs.BSplineCurve |
  bluemath.geom.nurbs.BezierSurface |
  bluemath.geom.nurbs.BSplineSurface;

declare function plot(data:PlotData|PlotData[], spec?:PlotSpec|PlotSpec[]);
`;

// Setup Monaco editor
require.config({ paths: { 'vs': './ext/monaco-editor/min/vs' }});
require(['vs/editor/editor.main'], function() {
  editor = monaco.editor.create(editorElem, {
    language: 'typescript',
    lineNumbers : false,
    minimap : {
      enabled : false
    }
  });

  // validation settings
  monaco.languages.typescript.typescriptDefaults.setDiagnosticsOptions({
    noSemanticValidation: true,
    noSyntaxValidation: false
  });

  // compiler options
  monaco.languages.typescript.typescriptDefaults.setCompilerOptions({
    target: monaco.languages.typescript.ScriptTarget.ES6,
    allowNonTsExtensions: true
  });

  // extra libraries
  monaco.languages.typescript.typescriptDefaults.addExtraLib(
    EXTRA_LIBS,'bluemath.d.ts');

  monaco.languages.typescript.typescriptDefaults.addExtraLib(
    preambleTypeDeclaration);

  // Setup example choices
  let exmNames = Object.keys(BMSHELL_EXAMPLES);
  let exmMap = {};
  for(let name of exmNames) {
    let key = nameToKey(name);
    let opt = document.createElement('option');
    opt.textContent = name;
    opt.setAttribute('value',key);
    exmChoices.appendChild(opt);
    exmMap[key] = name;
  }
  exmChoices.onchange = () => {
    let exmCode = BMSHELL_EXAMPLES[exmMap[exmChoices.value]];
    editor.getModel().setValue(exmCode);

    window.location.href = window.location.protocol + '//' +
      window.location.host + window.location.pathname + '#' +
      exmChoices.value;
  }

  // Populate with first example
  let urlmatch = /#([\d\w-]+)$/.exec(window.location.href);
  let selkey;
  if(urlmatch && exmMap.hasOwnProperty(urlmatch[1])) {
    selkey = urlmatch[1];
  } else {
    selkey = nameToKey(exmNames[0]);
  }
  exmChoices.value = selkey;
  editor.getModel().setValue(BMSHELL_EXAMPLES[exmMap[selkey]]);

  // Create IFrame
  ifrm = document.createElement('iframe');
  ifrm.setAttribute('src', './runtime.html');
  runtimeElem.appendChild(ifrm);
  ifrm.style.width = (W-500)+'px';
  ifrm.contentWindow.onload = function () {
    var runbutton = document.querySelector('#btn-run');
    runbutton.classList.remove('btn-disabled');
  }
});


// Recreate the runtime iframe and execute the code each time 'Run'
// is clicked
var runbutton = document.querySelector('#btn-run');
runbutton.classList.add('btn-disabled');
runbutton.onclick = () => {
  if(runbutton.classList.contains('btn-disabled')) {
    return; // done nothing
  }
  let codestring = editor.getModel().getValue();
  let ifrm = runtimeElem.querySelector('iframe');
  ifrm.contentWindow.scopedEval(codestring);
};

let docsButton = document.querySelector('#btn-docs');
docsButton.onclick = () => {
  window.open('http://www.bluemathsoftware.com/docs.html','_blank');
};

clearElem.onclick = () => {
  if(confirm('Clear the editor? You will loose any unsaved code')) {
    editor.getModel().setValue('');
  }
};

saveElem.onclick = () => {
  var blob = new Blob([editor.getModel().getValue()], {type: "text/plain;charset=utf-8"});
  saveAs(blob, 'bluemath-script.ts');
};

containerElem.style.width = W+'px';
containerElem.style.height = Math.round(0.9*H)+'px';
editorElem.style.width = '500px';
