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

const simple =
`
console.log('BlueMath version:',bluemath.version);
let A = new bluemath.common.NDArray([
  [1,0],
  [2,1],
  [3,4],
  [4,4],
  [5,8],
  [6,1]
]);
log('A',A.toString());
plot(A);
`;

const identity_and_multiplication =
`
const bm = bluemath;
let {eye,mul} = bm.common;

let A = eye(3,'i32');
log('Identity',A);
let B = mul(A,5);
let C = mul(A,2);
log('B',B);
log('C',C);
`;

const determinant = 
`
const {NDArray} = bluemath.common;
const {det,slogdet} = bluemath.linalg;

let A = new NDArray([
    [3,4,5],
    [0,3,4],
    [1,3,5]
]);

log('A',A);
log('Determinant = ',det(A));
let [sign,logdet] = slogdet(A);
log('Signed Log of Determinant = ',sign*logdet);
`;

const triangular_matrices = 
`
const {NDArray} = bluemath.common;
const {triu,tril} = bluemath.linalg;

let A = new NDArray([
    [4,5,6,7],
    [2,3,4,5],
    [0,2,3,4],
    [1,2,4,5]
],{datatype:'i32'});

plot(A,{type:'matrix',title:'A'});
plot([triu(A),tril(A)],
  [
    {type:'matrix',title:'Upper'},
    {type:'matrix',title:'Lower'}
  ]
);
`

const ndarray_arithmatic =
`
const {NDArray,eye,add,sub,mul,div} = bluemath.common;

let A = new NDArray([
    [2,3],
    [1,9]
],{datatype:'i32'});
let B = eye(2,'i32');
log('A',A);
log('B',B);
log('A+B',add(A,B));
log('A-B',sub(A,B));
log('3*A',mul(3,A));
log('A/2',div(A,2));
`

const ndarray_slicing =
`
const {NDArray} = bluemath.common;

let A = new NDArray([
    [1,2,3,4],
    [5,6,7,8],
    [9,10,11,12],
    [13,14,15,16]
],{datatype:'i32'});
log('A',A);
log('A[:,1] = ',A.get(':',1));
log('A[1,:] = ',A.get(1,':'));
log('A[null,2] = ',A.get(null,2));
log('A[:2,1:3] = ',A.get(':2','1:3'));
`

const nurbs =
`
let {arr} = bluemath.common;
let {
  BezierCurve,BSplineCurve,BezierSurface,BSplineSurface
} = bluemath.geom.nurbs;

let bezcrv = new BezierCurve(
    2,arr([[10,10],[30,25],[10,50]]));

let bcrv = new BSplineCurve(
    2,
    arr([
        [-5, -5],
        [-2, 0],
        [-1, 5],
        [-0.5, 2],
        [0.5, 2],
        [1, 5],
        [2, 0],
        [5, -5]
      ]),
    arr([0, 0, 0, 0.2, 0.4, 0.6, 0.8, 0.8, 1, 1, 1])
);

plot([bezcrv,bcrv],
    [{title:'Bezier Curve'}, {title:'BSpline Curve'}]);

let bezsrf = new BezierSurface(
    3, 2,
    [
      [[-1,-1,0],[0,-1,0],[1,-1,0]],
      [[-1,0,1],[0,0,2],[1,0,-1]],
      [[-1,2,0],[0,2,0],[1,2,0]],
      [[-1,4,0],[0,4,0],[1,4,0]]
    ]
);

let bsrf = new BSplineSurface(
  3, 2,
  [0,0,0,0,0.5,0.75,1,1,1,1],
  [0,0,0,0.5,1,1,1],
  [
    [ [-1,-1,2],[0,-1,1],[1,-1,1],[2,-1,1] ],
    [ [-1,0,1],[0,0,1],[1,0,1],[2,0,1] ],
    [ [-1,1,1],[0,1,-1],[1,1,-1],[2,1,-1] ],
    [ [-1,2,0],[0,2,0],[1,2,-1],[2,2,-1] ],
    [ [-1,3,0],[0,3,0],[1,3,-1],[2,3,-1] ],
    [ [-1,5,0],[0,5,0],[1,5,-1],[2,5,-1] ]
  ]
);
plot([bezsrf,bsrf],[
  {title:'Bezier Surface'},
  {title:'BSpline Surface'}
]);
`

let BMSHELL_EXAMPLES = {
  'Simple' : simple,
  'NDArray Arithmatic' : ndarray_arithmatic,
  'NDArray Slicing' : ndarray_slicing,
  'Identity and Multiplication' : identity_and_multiplication,
  'Matrix Determinant' : determinant,
  'Triangular Matrices' : triangular_matrices,
  'NURBS' : nurbs 
};
