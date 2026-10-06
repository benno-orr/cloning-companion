"""Manual macOS integration check: run with .venv/bin/python tests/native_drop_smoke.py.

Uses a separate test window; does not touch the user's running app. Simulates
the browser drop and native path handoff, not a physical Finder gesture.
"""
import json
from pathlib import Path
import sys
import tempfile
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import webview
from webview.dom import _dnd_state
from mac_app.main import NativeAPI, resource_path
from plasmid_verify.golden_gate_design import _write_snapgene_map


def main():
    api = NativeAPI()
    window = webview.create_window('CloningCompanion drop test', resource_path('mac_app/index.html').as_uri(), js_api=api)
    api.bind(window)
    errors = []

    def check():
        try:
            window.events.loaded.wait(15)
            deadline = time.monotonic() + 15
            while time.monotonic() < deadline:
                if window.evaluate_js('Boolean(window.appVersion)'):
                    break
                time.sleep(.1)
            assert window.evaluate_js('Boolean(window.appVersion)'), 'Window readiness never initialized app'
            assert _dnd_state['num_listeners'] > 0, 'No native drop listeners'
            assert api.register_design_drop_target(), 'Repeated registration failed'
            with tempfile.TemporaryDirectory(prefix='cloning-drop-test-') as directory:
                source = Path(directory) / 'annotated test.dna'
                _write_snapgene_map(source, 'CACC' + 'A' * 20 + 'TGAAATGGCC', [
                    {'name': '[Backbone]', 'start': 0, 'end': 28},
                    {'name': 'pUC ori', 'start': 4, 'end': 20},
                    {'name': '{Insert}', 'start': 28, 'end': 34},
                    {'name': '-', 'start': 24, 'end': 28},
                    {'name': '-', 'start': 0, 'end': 4},
                ], 'Drop integration test')
                # Cocoa normally obtains this mapping from the Finder pasteboard.
                _dnd_state['paths'].append((source.name, str(source)))
                window.evaluate_js('''(() => {
                    showView('design');
                    const transfer = new DataTransfer();
                    transfer.items.add(new File(['test'], %s));
                    document.getElementById('pick-design-map').dispatchEvent(
                        new DragEvent('drop', {bubbles:true, cancelable:true, dataTransfer:transfer}));
                })()''' % json.dumps(source.name))
                deadline = time.monotonic() + 15
                while time.monotonic() < deadline:
                    if window.evaluate_js('state.designMap?.path || null'):
                        break
                    time.sleep(.1)
                actual = window.evaluate_js('state.designMap')
                assert actual == api.annotated_design_map_from_path(str(source)), actual
                assert window.evaluate_js("document.getElementById('pick-design-map').classList.contains('loaded')")
                window.evaluate_js("showView('setup')")
                assert window.evaluate_js("document.body.dataset.workspace") == 'verification'
                assert window.evaluate_js("document.getElementById('design-navigation').classList.contains('hidden')")
                assert not window.evaluate_js("document.getElementById('verification-method').classList.contains('hidden')")
                window.evaluate_js("showView('design-results')")
                assert window.evaluate_js("document.body.dataset.workspace") == 'design'
                assert window.evaluate_js("document.getElementById('verification-method').classList.contains('hidden')")
                assert window.evaluate_js("document.getElementById('run-top').classList.contains('hidden')")
                assert window.evaluate_js('state.designMap') == actual, 'Switching workspaces discarded map'
                window.evaluate_js("document.getElementById('back-to-design').click()")
                assert window.evaluate_js('activeView') == 'design'
                assert not window.evaluate_js("document.getElementById('run-design-top').disabled")
                print('PASS: real WebKit readiness → native listener → browser drop → path resolution → map import → loaded UI', flush=True)
                print('PASS: separate workspace navigation, contextual actions, empty outputs, and retained design inputs', flush=True)
                window.evaluate_js("document.getElementById('check-updates').click()")
                deadline = time.monotonic() + 10
                while time.monotonic() < deadline:
                    if 'Installed:' in window.evaluate_js("document.getElementById('modal-body').innerText"):
                        break
                    time.sleep(.1)
                assert 'Installed:' in window.evaluate_js("document.getElementById('modal-body').innerText")
                assert not window.evaluate_js("document.getElementById('modal-backdrop').classList.contains('hidden')")
                window.evaluate_js("document.getElementById('modal-close').click()")
                assert window.evaluate_js("document.getElementById('modal-backdrop').classList.contains('hidden')")
                print('PASS: in-app update button → native check → version/status dialog', flush=True)
                designed = api.run_annotated_golden_gate_design(str(source), str(Path(directory)/'outputs'))
                assert designed['ok'] and designed['plasmidCount'] == 1
                window.evaluate_js('renderDesignResults(%s)' % json.dumps(designed))
                assert window.evaluate_js("document.querySelectorAll('[data-plasmid]').length") == 1
                assert 'Assembled plasmids' in window.evaluate_js("document.getElementById('design-results').innerText")
                print('PASS: unique SnapGene outputs appear in app with per-file actions', flush=True)
        except Exception as exc:
            errors.append(exc)
            print('FAIL:', repr(exc), flush=True)
        finally:
            window.destroy()

    webview.start(check, private_mode=True)
    if errors:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
