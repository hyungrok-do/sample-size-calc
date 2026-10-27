import { mountDag } from './dag.js?v=dev';

mountDag(document.getElementById('dag'));
if ('serviceWorker' in navigator) navigator.serviceWorker.register('sw.js').catch(() => {});
