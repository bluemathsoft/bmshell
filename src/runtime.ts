
import * as bluemath from 'bluemath'

(<any>window).bluemath = bluemath;

(<any>window).console = {
  log : function(...args:any[]) {
    let msg:string = args.map(a => new String(a)).join('');
    let el = document.createElement('p');
    el.textContent = msg;
    document.body.appendChild(el);
  }
};