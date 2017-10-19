
const fs = require('fs');
const path = require('path');

let node_modules = path.join(__dirname,'..','node_modules');

let outfile = fs.openSync('build/extralibs.js','w');

fs.writeSync(outfile,'var EXTRA_LIBS = `\n');

let dirstack = [];
let namespace;

function walk(topdir) {
  let dirpath = path.join(topdir,...dirstack);
  let fnames = fs.readdirSync(dirpath);
  for(let fname of fnames) {
    let stat = fs.statSync(path.join(dirpath,fname));
    if(stat.isDirectory()) {
      dirstack.push(fname);
      namespace.push(fname);
      walk(topdir);
      namespace.pop();
      dirstack.pop();
    } else if(fname.endsWith('.d.ts')) {
      let fpath = path.join(dirpath,fname);
      let source = new String(fs.readFileSync(fpath));
      for(let line of source.split('\n')) {
        if(/import/g.test(line)) {
        } else if(/`/g.test(line)) {
        } else if(/export\s+(?!(declare))/g.test(line)) {
        } else {
          line = line.replace(/export declare/g,'declare');
          fs.writeSync(outfile,'    '+line+'\n');
        }
      }
    }
  }
}

namespace = ['bluemath']

fs.writeSync(outfile,'declare namespace bluemath {\n');

fs.writeSync(outfile,'  declare namespace common {\n');
walk(path.join(node_modules,'@bluemath','common','lib'));
fs.writeSync(outfile,'  }\n'); // namespace common 

fs.writeSync(outfile,'  declare namespace linalg {\n');
walk(path.join(node_modules,'@bluemath','linalg','lib'));
namespace.pop();
fs.writeSync(outfile,'  }\n'); // namespace linalg 

fs.writeSync(outfile,'  declare namespace geom {\n');
walk(path.join(node_modules,'@bluemath','geom','lib'));
fs.writeSync(outfile,'  }\n'); // namespace linalg 

fs.writeSync(outfile,'  declare namespace topo {\n');
walk(path.join(node_modules,'@bluemath','topo','lib'));
fs.writeSync(outfile,'  }\n'); // namespace linalg 

fs.writeSync(outfile,'}\n'); // namespace bluemath

fs.writeSync(outfile,'`;\n');