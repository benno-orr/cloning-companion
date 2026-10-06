// Regression: the desktop framework dispatches readiness on window, not document.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const html = fs.readFileSync('mac_app/index.html', 'utf8');
for (const match of html.matchAll(/<script[^>]*>([\s\S]*?)<\/script>/g)) new vm.Script(match[1]);
const bootstrap = html.slice(html.indexOf('let nativeBridgeInitialization='), html.indexOf('// Initial workspace render.'));

async function check(alreadyReady) {
  const window = new EventTarget();
  let registrations = 0, renders = 0;
  const api = {
    register_design_drop_target: async () => { registrations++; return true; },
    app_info: async () => ({version: 'test', recent: []}),
  };
  const context = vm.createContext({window, state: {}, renderRecent() {renders++;}, renderInputs() {},
    toast(message) {throw new Error(message);}});
  if (alreadyReady) window.pywebview = context.pywebview = {api};
  vm.runInContext(bootstrap, context);
  if (!alreadyReady) {
    assert.equal(registrations, 0);
    window.pywebview = context.pywebview = {api};
  }
  window.dispatchEvent(new Event('pywebviewready'));
  window.dispatchEvent(new Event('pywebviewready'));
  await vm.runInContext('initializeNativeBridge()', context);
  assert.equal(registrations, 1, 'drop listener registered once');
  assert.equal(renders, 1);
  assert.equal(window.appVersion, 'test');
}
(async () => {
  await check(false);
  await check(true);
  console.log('Native bridge readiness tests passed (window event, early bridge, duplicate events).');
})().catch(error => {console.error(error); process.exitCode = 1;});
