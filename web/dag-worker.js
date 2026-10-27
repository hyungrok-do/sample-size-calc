import { simulateGraph } from './dag-model.js?v=dev';

self.onmessage = ({ data }) => {
  try {
    const result = simulateGraph(data.graph, data.options);
    self.postMessage({ result }, Object.values(result.columns).map((column) => column.buffer));
  } catch (error) {
    self.postMessage({ error: error.message });
  }
};
