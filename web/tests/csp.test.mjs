// The page's Content-Security-Policy in the source: the <style>'s hash, and the same policy in the service worker.
import assert from 'node:assert/strict';
import test from 'node:test';
import { readFile } from 'node:fs/promises';
import { policy } from './csp.mjs';

const read = (file) => readFile(new URL(`../${file}`, import.meta.url), 'utf8');
const [html, sw] = await Promise.all([read('index.html'), read('sw.js')]);

test('the policy lets in the page\'s <style> by its hash, and sw.js sends the same policy', () => {
  const { meta } = policy(html, sw);
  assert.match(meta, /^default-src 'none';/);
});

test('an edit to the <style> without a new hash is caught', () => {
  assert.throws(() => policy(html.replace('--bg: #f6f6f3;', '--bg: #f6f6f4;'), sw), /hashes to/);
  assert.throws(() => policy(html.replace('<style>\n', '<style>\n\n'), sw), /hashes to/);
  assert.doesNotThrow(() => policy(html.replaceAll('\n', '\r\n'), sw));  // as a browser reads a file saved with CRLF
});

test('a policy in sw.js that differs from the page\'s is caught', () => {
  assert.throws(() => policy(html, sw.replace("form-action 'none'; frame-ancestors", 'frame-ancestors')), /not the page's/);
  assert.throws(() => policy(html, sw.replace("headers.set('Content-Security-Policy', POLICY);", '')), /as a header/);
  assert.throws(() => policy(html.replace('<button id="skip"', '<button style="color: red" id="skip"'), sw), /blocks/);
});
