// This file houses the code for all the interactive data and charts on the dashboard.

//--------------------------------------------
// Helper functions
//--------------------------------------------

// Function that copies uncorrectable barcodes to clipboard
function copyUncorrectableBarcodes(type) {
    if (!data.top_uncorrectables[type]) {
        return alert("No data available for " + type);
    }
    navigator.clipboard.writeText(data.top_uncorrectables[type].map(d => `${d.barcode}`).join("\n"));
}

// Attach click listeners to all .copy-btn buttons
document.querySelectorAll(".copy-btn").forEach(btn => {
    btn.addEventListener("click", () => {
        copyUncorrectableBarcodes(btn.dataset.type);
    });
});

// Hide charts for ignored barcodes.
document.addEventListener("DOMContentLoaded", function () {
  if (data.ignore_p5) {
    document.getElementById("container-well-p5").style.display = "none";
    document.getElementById("container-uncorrectables-p5").style.display = "none";
  }
  if (data.ignore_p7) {
    document.getElementById("container-well-p7").style.display = "none";
    document.getElementById("container-uncorrectables-p7").style.display = "none";
  }
});

// Helper function to create an element with class and innerHTML
const createElement = (type, className, innerHTML) => {
  const element = document.createElement(type);
  element.className = className;
  if (innerHTML) element.innerHTML = innerHTML;
  return element;
  };

function initializePopovers() {
  const warningText = "Saturation below 40%: take cell recovery with a grain of salt";
  const warningElements = document.querySelectorAll(".saturation-warning");

  if (window.tippy) {
    warningElements.forEach((el) => {
      if (el._tippy) el._tippy.destroy();
      tippy(el, {
        content: warningText,
        trigger: "mouseenter focus",
        placement: "top",
      });
    });
    return;
  }

  warningElements.forEach((el) => {
    el.setAttribute("title", warningText);
  });
}


// Function to generate a matrix chart using heatmap.js
function generateMatrixChart(data_chart, ctx, rgb_x, rgb_y, rgb_z) {

    // Convert the data to a format that heatmap.js can read.
    var data_index = [];

    // Determine max frequency of data_chart
    var max_frequency = 0;
    for (var i = 0; i < data_chart.length; i++) {
      if (data_chart[i].frequency > max_frequency) {
        max_frequency = data_chart[i].frequency;
      }
    }

    for (var i = 0; i < data_chart.length; i++) {
      data_index.push({
        x: data_chart[i].col,
        y: data_chart[i].row,
        v: data_chart[i].frequency,
        v_perc: data_chart[i].frequency / max_frequency
      });
    }

    new Chart(ctx, {
      type: 'matrix',
      data: {
        datasets: [{
          data: data_index,
          // Color alpha based on frequency compared to max frequency
          backgroundColor(context) {
            var v_perc = context.dataset.data[context.dataIndex].v_perc;
            return 'rgba(' + rgb_x + ', ' + rgb_y + ', ' + rgb_z + ', ' + v_perc + ')';
          },
          borderColor() {
            return 'rgb(0, 0, 0)';
          },
          // Round the corners of the squares
          borderRadius: 2,
          borderWidth: .5,
          width(context) {
            const a = context.chart.chartArea;
            if (!a) {
              return 0;
            }
            return (a.right - a.left) / 12 - 2;
          },
          height(context) {
            const a = context.chart.chartArea;
            if (!a) {
              return 0;
            }
            return (a.bottom - a.top) / 8 - 2;
          }
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          legend: false,
          tooltip: {
            callbacks: {
              title() {
                return '';
              },
              label(context) {
                const v = context.dataset.data[context.dataIndex];
                return [Intl.NumberFormat("en-US").format(v.v) + ' (' + v.y + v.x + ')'];
              }
            }
          }
        },
        scales: {
          x: {
            type: 'category',
            labels: ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'],
            offset: true,
            ticks: {
              display: true
            },
            grid: {
              display: false
            },
          },
          y: {
            type: 'category',
            labels: ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
            offset: true,
            reverse: false,
            ticks: {
              display: true
            },
            grid: {
              display: false
            },
          }
        }
      }
    }
    );
}

//--------------------------------------------
// Insert data from qc_data.js
//--------------------------------------------

const sample_names = Object.keys(data.sample_success);
const roundToOne = num => +(Math.round(num + "e+1") + "e-1");

// Count the total number of estimated cells over samples.
let total_estimated_cells = 0;
for (let i = 0; i < sample_names.length; i++) {
  total_estimated_cells += data.sample_success[sample_names[i]].estimated_cells;
}

// Update numbers
document.addEventListener("DOMContentLoaded", () => {
  document.getElementById("version").innerHTML = data.version;
  document.getElementById("experiment_name").innerHTML = data.experiment_name;
  document.getElementById("n_totalsamples").innerHTML = Intl.NumberFormat("en-US").format(sample_names.length);
  document.getElementById("n_total_pairs").innerHTML = Intl.NumberFormat("en-US").format(data.n_pairs);
  document.getElementById("n_total_pairs_success_perc").innerHTML = `${roundToOne((data.n_pairs_success / data.n_pairs) * 100)}%`;
  document.getElementById("n_total_pairs_failure_perc").innerHTML = `${roundToOne((data.n_pairs_failure / data.n_pairs) * 100)}%`;
  document.getElementById("n_total_corrections").innerHTML = Intl.NumberFormat("en-US").format((data.ignore_p5 ? 0 : data.n_corrected_p5) + (data.ignore_p7 ? 0 : data.n_corrected_p7) + data.n_corrected_ligation + data.n_corrected_rt + data.n_corrected_hashing);
  document.getElementById("n_total_cells").innerHTML = Intl.NumberFormat("en-US").format(total_estimated_cells);

  // Set the n_total_pairs_success_perc_bar width
  document.getElementById("n_total_pairs_success_perc_bar").style.width = `${(data.n_pairs_success / data.n_pairs) * 100}%`;
});

//--------------------------------------------
// Chart - Rescues overview (pie-chart)
//--------------------------------------------

// Generate the data for the doughnut chart.
function buildCorrectableBarcodesDoughnutData(data) {
  const specs = [
    { label: "p5", key: "n_corrected_p5", color: "#d63939", ignore: data.ignore_p5 },
    { label: "p7", key: "n_corrected_p7", color: "#1f77b4", ignore: data.ignore_p7 },
    { label: "ligation", key: "n_corrected_ligation", color: "#ff7f0e" },
    { label: "RT", key: "n_corrected_rt", color: "#2ca02c" },
    { label: "hashing", key: "n_corrected_hashing", color: "#9467bd" },
  ];

  const labels = [];
  const values = [];
  const colors = [];

  for (const { label, key, color, ignore } of specs) {
    // Skip ignored barcodes or zero/invalid counts.
    if (ignore) continue;
    const v = Number(data?.[key] ?? 0);
    if (!Number.isFinite(v) || v === 0) continue;

    labels.push(label);
    values.push(v);
    colors.push(color);
  }

  return {
    labels,
    datasets: [
      {
        data: values,
        backgroundColor: colors,
      },
    ],
  };
}

// Generate the doughnut chart.
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-rescues").innerHTML = '';
  var ctx = document.getElementById("chart-rescues").getContext("2d");
  var chart = new Chart(ctx, {
    type: "doughnut",
    data: buildCorrectableBarcodesDoughnutData(data),
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: {
          position: "right",
          labels: {
            boxWidth: 15,
            font: {
              size: 10,
            },
          },
        },
      },
    },
  });
});

//--------------------------------------------
// Chart - No. of successful pairs per sample.
//--------------------------------------------

var sample_n_pairs_success = [];
for (var i = 0; i < sample_names.length; i++) {
  sample_n_pairs_success.push({
    sample_name: sample_names[i],
    frequency: data.sample_success[sample_names[i]].n_pairs_success,
  });
}

// Sort the array by frequency
sample_n_pairs_success.sort(function (a, b) {
  return b.frequency - a.frequency;
});


document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-n_pairs_sample").innerHTML = ''
  var ctx = document.getElementById("chart-n_pairs_sample").getContext("2d");

  const chart_n_pair = new Chart(ctx, {
    type: "bar",
    data: {
      labels: sample_n_pairs_success.map(function (d) {
        return d.sample_name;
      }),
      datasets: [
        {
          label: "No. of successful pairs",
          data: sample_n_pairs_success.map(function (d) {
            return d.frequency;
          }),
          backgroundColor: [
            'rgba(255, 99, 132, 0.9)',
            'rgba(255, 159, 64, 0.9)',
            'rgba(255, 205, 86, 0.9)',
            'rgba(75, 192, 192, 0.9)',
            'rgba(54, 162, 235, 0.9)',
            'rgba(153, 102, 255, 0.9)',
            'rgba(201, 203, 207, 0.9)'
          ],
          borderColor: [
            'rgba(0, 0, 0)',
          ],
          borderWidth: 1,
        },
      ],
    },
    options: {
      indexAxis: 'x',
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: {
          display: false,
        },
      },
      scales: {
        x: {
          grid: {
            display: false
          },
          ticks: {
            maxRotation: 90,
            minRotation: 90,
          }
        },
      }
    },
  });
});

//--------------------------------------------
// Heatmaps - Index counts
//--------------------------------------------

// Well - p5
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-well-p5").innerHTML = ''
  var ctx = document.getElementById("chart-well-p5").getContext("2d");
  generateMatrixChart(data.p5_index_counts, ctx, 31, 119, 180);
});

// Well - p7
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-well-p7").innerHTML = ''
  var ctx = document.getElementById("chart-well-p7").getContext("2d");
  generateMatrixChart(data.p7_index_counts, ctx, 255, 105, 180);
});

// Well - RT plate 1
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-rt_plate_01").innerHTML = ''
  var ctx = document.getElementById("chart-rt_plate_01").getContext("2d");

  // Display "No data" if there if undefined.
  if (data.rt_barcode_counts.P01 == undefined) {
    document.getElementById("rt_plate_01").innerHTML = "<div style='text-align:center'><b>Plate not used</b></div>";
  }else{
    generateMatrixChart(data.rt_barcode_counts.P01, ctx, 255, 140, 0);
  }
});

// Well - RT plate 2
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-rt_plate_02").innerHTML = ''
  var ctx = document.getElementById("chart-rt_plate_02").getContext("2d");

  // Display "No data" if there if undefined.
  if (data.rt_barcode_counts.P02 == undefined) {
    document.getElementById("rt_plate_02").innerHTML = "<div style='text-align:center'><b>Plate not used</b></div>";

  }else{
    generateMatrixChart(data.rt_barcode_counts.P02, ctx, 255, 140, 0);
  }
});

// Well - RT plate 3
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-rt_plate_03").innerHTML = ''
  var ctx = document.getElementById("chart-rt_plate_03").getContext("2d");
    // Display "No data" if there if undefined.
  if (data.rt_barcode_counts.P03 == undefined) {
    document.getElementById("rt_plate_03").innerHTML = "<div style='text-align:center'><b>Plate not used</b></div>";

  }else{
    generateMatrixChart(data.rt_barcode_counts.P03, ctx, 255, 140, 0);
  }
});

// Well - RT plate 4
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-rt_plate_04").innerHTML = ''
  var ctx = document.getElementById("chart-rt_plate_04").getContext("2d");

  if (data.rt_barcode_counts.P04 == undefined) {
    document.getElementById("rt_plate_04").innerHTML = "<div style='text-align:center'><b>Plate not used</b></div>";
  }else{
    generateMatrixChart(data.rt_barcode_counts.P04, ctx, 255, 140, 0);
  }

});

//--------------------------------------------
// Chart - Ligation usage.
//--------------------------------------------

document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-ligation").innerHTML = ''
  var ctx = document.getElementById("chart-ligation").getContext("2d");

  // Order on name LIG01, LIG02, LIG03, LIG04, LIG05, LIG06, LIG07, LIG08, LIG09, LIG10, LIG11, LIG12, etc.
  var ligation_barcode_counts = [];
  for (var i = 0; i < data.ligation_barcode_counts.length; i++) {
    ligation_barcode_counts.push({
      barcode: data.ligation_barcode_counts[i].barcode,
      frequency: data.ligation_barcode_counts[i].frequency,
      // Remove 0 from the index if it is present.
      index: data.ligation_barcode_counts[i].barcode.replace("LIG", "").replace("LIG0", ""),
    });
  }

  // Sort the array by index
  ligation_barcode_counts.sort(function (a, b) {
    return a.index - b.index;
  });


  const chart_ligation = new Chart(ctx, {
    type: "bar",
    data: {
      labels: ligation_barcode_counts.map(function (d) {
        return d.barcode;
      }),
      datasets: [
        {
          label: "Frequency",
          data: ligation_barcode_counts.map(function (d) {
            return d.frequency;
          }),
          // Color the first 96 barcodes in red, the next 96 in blue, the next 96 in green and the last 96 in yellowu sing context.
          backgroundColor(context) {
            if (context.dataIndex < 96) {
              return 'rgba(255, 99, 132, 0.9)';
            } else if (context.dataIndex < 192) {
              return 'rgba(54, 162, 235, 0.9)';
            } else if (context.dataIndex < 288) {
              return 'rgba(255, 205, 86, 0.9)';
            } else {
              return 'rgba(75, 192, 192, 0.9)';
            }
          },
          borderColor: [
            'rgba(0, 0, 0)',
          ],
          borderWidth: .1,
        },
      ],
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: {
          display: false,
        },
      },
      scales: {
        x: {
          grid: {
            display: false
          },
          ticks: {
            maxRotation: 90,
            minRotation: 90,
          }
        },
      }
    },
  });
});

//--------------------------------------------
// Chart - Top 15 uncorrectable barcodes.
//--------------------------------------------


function generateChart_uncorrectables(id, data, color) {
  document.addEventListener("DOMContentLoaded", function () {
    document.getElementById(id).innerHTML = '';
    var ctx = document.getElementById(id).getContext("2d");

    // Order on frequency
    var chartData = [];
    for (var i = 0; i < data.length; i++) {
      chartData.push({
        barcode: data[i].barcode,
        frequency: data[i].frequency,
      });
    }

    // Sort the array by frequency
    chartData.sort(function (a, b) {
      return b.frequency - a.frequency;
    });

    new Chart(ctx, {
      type: "bar",
      data: {
        labels: chartData.map(function (d) {
          return d.barcode;
        }),
        datasets: [
          {
            label: "Frequency",
            data: chartData.map(function (d) {
              return d.frequency;
            }),
            backgroundColor: color,
            borderColor: [
              'rgba(0, 0, 0)',
            ],
            borderWidth: 1,
          },
        ],
      },
      options: {
        indexAxis: 'y',
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          legend: {
            display: false,
          },
        },
        scales: {
          y: {
            ticks: {
              font: {
                family: "Courier New",
                size: 12,
                autoSkip: false,
              },
            },
          },
          x: {
            grid: {
              display: false
            },
            ticks: {
              font: {
                family: "Courier New",
                size: 10,
              },
            },
          }
        }
      },
    });
  });
}

// Usage example:
generateChart_uncorrectables("chart-top_uncorrectables_p5", data.top_uncorrectables.p5, color = "#d63939");
generateChart_uncorrectables("chart-top_uncorrectables_p7", data.top_uncorrectables.p7, color = "#1f77b4");
generateChart_uncorrectables("chart-top_uncorrectables_lig", data.top_uncorrectables.ligation, color = "#ff7f0e");
generateChart_uncorrectables("chart-top_uncorrectables_rt", data.top_uncorrectables.rt, color = "#2ca02c");

//--------------------------------------------
// Export data to JSON table.
//--------------------------------------------

// Print data to json-table (div)
document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("json-data").innerHTML = JSON.stringify(data, undefined, 2);
});

//--------------------------------------------
// Table - STARSolo summary.
//--------------------------------------------

function generate_starsolo_table(data) {
  const table = document.getElementById("sample-starsolo-table");
  table.className = "table card-table table-vcenter text-nowrap datatable tablesorter";
  table.innerHTML = `
      <thead>
        <tr>
          <th>Sample</th>
          <th>Input reads</th>
          <th>Sequencing saturation</th>
          <th>Mapped reads<br><sub>(Genome; %)</sub></th>
          <th>Mapped reads<br><sub>(Genome - Unique; %)</sub></th>
          <th>Mapped reads<br><sub>(Genes; %)</sub></th>
          <th>Mapped reads<br><sub>(Genes - Unique; %)</sub></th>
          <th># intronic reads</th>
          <th># intergenic reads</th>
          <th># mitochondrial reads</th>
          <th># exonic (AS) reads</th>
          <th># intronic (AS) reads</th>
          <th>Reads per cell<br><sub>(Mean)</sub></th>
          <th>Genes per cell<br><sub>(Mean)</sub></th>
          <th>UMI per cell<br><sub>(Mean)</sub></th>
          <th>Estimated cells</th>
        </tr>
      </thead>
      <tbody class="table-tbody">
      </tbody>`;

  const tbody = table.querySelector("tbody");
  for (const sample in data) {
    const sequencingSaturationPercent = Math.round(data[sample].sequencing_saturation * 100);
    let sequencingSaturationTextClass = "";
    let sequencingSaturationBarClass = "";
    if (sequencingSaturationPercent < 20) {
      sequencingSaturationTextClass = "text-red";
      sequencingSaturationBarClass = "bg-red";
    } else if (sequencingSaturationPercent < 40) {
      sequencingSaturationTextClass = "text-yellow";
      sequencingSaturationBarClass = "bg-yellow";
    }

    const row = tbody.insertRow(-1);
    row.insertCell(0).innerHTML = sample;
    row.insertCell(1).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_reads);
    const saturationCell = row.insertCell(2);
    const progress = createElement("div", "row align-items-center");
    const col1 = createElement("div", `col-12 col-lg-auto ${sequencingSaturationTextClass}`, `${sequencingSaturationPercent}%`);
    if (sequencingSaturationPercent < 40) {
      col1.classList.add("saturation-warning");
      col1.setAttribute("tabindex", "0");
      col1.style.cursor = "help";
    }
    const col2 = createElement("div", "col");
    const progress_bar = createElement("div", "progress");
    progress_bar.style.width = "3rem";
    const progress_bar_inner = createElement("div", `progress-bar ${sequencingSaturationBarClass}`);
    progress_bar.style.height = "0.55rem";
    if (sequencingSaturationPercent < 40) {
      progress_bar.classList.add("saturation-warning");
      progress_bar.setAttribute("tabindex", "0");
      progress_bar.style.cursor = "help";
      progress_bar_inner.classList.add("saturation-warning");
      progress_bar_inner.setAttribute("tabindex", "0");
      progress_bar_inner.style.cursor = "help";
    }
    progress_bar_inner.style.width = `${sequencingSaturationPercent}%`;
    progress_bar_inner.setAttribute("role", "progressbar");
    progress_bar_inner.setAttribute("aria-valuenow", sequencingSaturationPercent);
    progress_bar_inner.setAttribute("aria-valuemin", "0");
    progress_bar_inner.setAttribute("aria-valuemax", "100");
    const span = createElement("span", "visually-hidden", `${sequencingSaturationPercent}%`);
    progress_bar_inner.appendChild(span);
    progress_bar.appendChild(progress_bar_inner);
    col2.appendChild(progress_bar);
    progress.appendChild(col1);
    progress.appendChild(col2);
    saturationCell.appendChild(progress);

    row.insertCell(3).innerHTML = `${roundToOne(data[sample].perc_mapped_reads_genome * 100)}%`;
    row.insertCell(4).innerHTML = `${roundToOne(data[sample].perc_unique_reads_genome_unique * 100)}%`;
    row.insertCell(5).innerHTML = `${roundToOne(data[sample].perc_mapped_reads_gene * 100)}%`;
    row.insertCell(6).innerHTML = `${roundToOne(data[sample].perc_unique_reads_gene_unique * 100)}%`;
    row.insertCell(7).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_intronic_reads);
    row.insertCell(8).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_intergenic_reads);
    row.insertCell(9).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_mitochondrial_reads);
    row.insertCell(10).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_exonicAS_reads);
    row.insertCell(11).innerHTML = Intl.NumberFormat("en-US").format(data[sample].total_intronicAS_reads);
    row.insertCell(12).innerHTML = Intl.NumberFormat("en-US").format(data[sample].mean_reads_per_cell);
    row.insertCell(13).innerHTML = Intl.NumberFormat("en-US").format(data[sample].mean_genes_per_cell);
    row.insertCell(14).innerHTML = Intl.NumberFormat("en-US").format(data[sample].mean_umi_per_cell);
    row.insertCell(15).innerHTML = Intl.NumberFormat("en-US").format(data[sample].estimated_cells);
  }
}

document.addEventListener("DOMContentLoaded", function () {
  generate_starsolo_table(data.sample_success);
  initializePopovers();

  // Sort using tablesorter
  $("#sample-starsolo-table").tablesorter({
    sortList: [[0, 0]],
  });

});

//--------------------------------------------
// Table - Hashing summary.
//--------------------------------------------

function insertFractionPassingCell(row, fractionPassing) {
  const fractionCell = row.insertCell(-1);
  if (fractionPassing == null) {
    fractionCell.innerHTML = "NA";
    return;
  }

  const fractionPassingPercent = Math.round(Number(fractionPassing) * 100);
  let fractionPassingTextClass = "";
  let fractionPassingBarClass = "";
  if (fractionPassingPercent < 20) {
    fractionPassingTextClass = "text-red";
    fractionPassingBarClass = "bg-red";
  } else if (fractionPassingPercent < 40) {
    fractionPassingTextClass = "text-yellow";
    fractionPassingBarClass = "bg-yellow";
  }

  const progress = createElement("div", "row align-items-center");
  const col1 = createElement("div", `col-12 col-lg-auto ${fractionPassingTextClass}`, `${fractionPassingPercent}%`);
  const col2 = createElement("div", "col");
  const progress_bar = createElement("div", "progress");
  progress_bar.style.width = "3rem";
  progress_bar.style.height = "0.55rem";
  const progress_bar_inner = createElement("div", `progress-bar ${fractionPassingBarClass}`);
  progress_bar_inner.style.width = `${fractionPassingPercent}%`;
  progress_bar_inner.setAttribute("role", "progressbar");
  progress_bar_inner.setAttribute("aria-valuenow", fractionPassingPercent);
  progress_bar_inner.setAttribute("aria-valuemin", "0");
  progress_bar_inner.setAttribute("aria-valuemax", "100");
  const span = createElement("span", "visually-hidden", `${fractionPassingPercent}%`);
  progress_bar_inner.appendChild(span);
  progress_bar.appendChild(progress_bar_inner);
  col2.appendChild(progress_bar);
  progress.appendChild(col1);
  progress.appendChild(col2);
  fractionCell.appendChild(progress);
}

function generate_hashing_summary_table(tableId, summaryRows) {
  const table = document.getElementById(tableId);
  if (!table) return;

  table.className = "table card-table table-vcenter text-nowrap datatable tablesorter";
  table.innerHTML = `
      <thead>
        <tr>
          <th>Sample</th>
          <th>Mean hash count</th>
          <th>Hash count total</th>
          <th>Hash UMI total</th>
          <th>Median ratio</th>
          <th>Mean ratio</th>
          <th>Total cells</th>
          <th>Cells passing<br><sub>(ratio >= 3)</sub></th>
          <th>Fraction passing</th>
        </tr>
      </thead>
      <tbody class="table-tbody">
      </tbody>`;

  const tbody = table.querySelector("tbody");
  const formatRatioValue = (value) => {
    if (value == null) return "NA";
    if (typeof value === "number" && !Number.isFinite(value)) {
      return value > 0 ? "Inf" : "-Inf";
    }
    return roundToOne(value);
  };

  for (const rowData of summaryRows) {
    const meanHashCount = rowData.mean_hash_count;
    const hashCountTotal = rowData.hash_count_total;
    const hashUmiTotal = rowData.hash_umi_total;
    const medianHashRatio = rowData.median_hash_ratio;
    const meanHashRatio = rowData.mean_hash_ratio;
    const totalCells = rowData.total_cells;

    const row = tbody.insertRow(-1);
    row.insertCell(0).innerHTML = rowData.sample_name ?? "";
    row.insertCell(1).innerHTML = meanHashCount == null ? "NA" : roundToOne(meanHashCount);
    row.insertCell(2).innerHTML = hashCountTotal == null ? "NA" : Intl.NumberFormat("en-US").format(hashCountTotal);
    row.insertCell(3).innerHTML = hashUmiTotal == null ? "NA" : Intl.NumberFormat("en-US").format(hashUmiTotal);
    row.insertCell(4).innerHTML = formatRatioValue(medianHashRatio);
    row.insertCell(5).innerHTML = formatRatioValue(meanHashRatio);
    row.insertCell(6).innerHTML = totalCells == null ? "NA" : Intl.NumberFormat("en-US").format(totalCells);
    row.insertCell(7).innerHTML = rowData.cells_passing == null ? "NA" : Intl.NumberFormat("en-US").format(rowData.cells_passing);
    insertFractionPassingCell(row, rowData.fraction_passing);
  }
}

function generate_hashing_bins_heatmap(binRows, binLabels) {
  const container = document.getElementById("chart-hashing-bins-heatmap");
  if (!container) return;

  if (typeof Plotly === "undefined") {
    container.innerHTML = "<div style='text-align:center'><b>Heatmap library unavailable</b></div>";
    return;
  }

  const xLabels = Array.isArray(binLabels) && binLabels.length > 0
    ? binLabels
    : [...new Set((binRows || []).map((row) => row.count_bin))];
  const yLabels = [...new Set((binRows || []).map((row) => row.sample_name))];

  if (xLabels.length === 0 || yLabels.length === 0) {
    container.innerHTML = "<div style='text-align:center'><b>No bin data available</b></div>";
    return;
  }

  const xIndexByLabel = Object.fromEntries(xLabels.map((label, idx) => [label, idx]));
  const yIndexByLabel = Object.fromEntries(yLabels.map((label, idx) => [label, idx]));
  const zMatrix = yLabels.map(() => xLabels.map(() => null));
  const nCellMatrix = yLabels.map(() => xLabels.map(() => 0));

  for (const row of (binRows || [])) {
    if (!row || !(row.count_bin in xIndexByLabel) || !(row.sample_name in yIndexByLabel)) continue;
    const xIdx = xIndexByLabel[row.count_bin];
    const yIdx = yIndexByLabel[row.sample_name];
    nCellMatrix[yIdx][xIdx] = Number(row.n_cells ?? 0);

    if (row.n_cells > 0 && row.fraction_passing != null) {
      zMatrix[yIdx][xIdx] = Number(row.fraction_passing) * 100;
    }
  }

  const hasData = zMatrix.some((row) => row.some((value) => value != null && Number.isFinite(value)));
  if (!hasData) {
    container.innerHTML = "<div style='text-align:center'><b>No bin data available</b></div>";
    return;
  }

  const trace = {
    type: "heatmap",
    x: xLabels,
    y: yLabels,
    z: zMatrix,
    customdata: nCellMatrix,
    colorscale: "Viridis",
    zmin: 0,
    zmax: 100,
    colorbar: {
      title: {
        text: "Passing (%)",
      },
      ticksuffix: "%",
    },
    hovertemplate: "%{y} | %{x}<br>Passing: %{z:.1f}%<br>n=%{customdata}<extra></extra>",
    hoverongaps: false,
  };

  const layout = {
    autosize: true,
    dragmode: false,
    hovermode: "closest",
    margin: { l: 120, r: 20, t: 20, b: 100 },
    xaxis: {
      tickangle: -45,
      automargin: true,
      title: { text: "Hash count bin" },
      fixedrange: true,
    },
    yaxis: {
      automargin: true,
      title: { text: "Sample" },
      fixedrange: true,
    },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(0,0,0,0)",
  };

  const config = {
    responsive: true,
    displayModeBar: false,
    scrollZoom: false,
    doubleClick: false,
    displaylogo: false,
  };

  Plotly.react(container, [trace], layout, config);
}

// For each hashing barcode, generate a row in the table.
// The row will contain the following information:
// 1. Hashing Sample
// 2. Hashing barcode
// 3. Total count (summed)
// 4. Total count (correct)
// 5. Total count (corrected)
// 6. Total count (correct - upstream)
function generate_hashing_table(data) {
  const table = document.getElementById("sample-hashing-table");
  table.className = "table card-table table-vcenter text-nowrap datatable tablesorter";
  table.innerHTML = `
      <thead>
        <tr>
          <th>Sample</th>
          <th>Hashing barcode</th>
          <th>Total count<br><sub>(summarized)</sub></th>
          <th>Total count<br><sub>(correct)</sub></th>
          <th>Total count<br><sub>(corrected)</sub></th>
          <th>Total count<br><sub>(correct:upstream)</sub></th>
        </tr>
      </thead>
      <tbody class="table-tbody">
      </tbody>`;
  const tbody = table.querySelector("tbody");
  for (var sample in data) {
    for (var hashing_barcode in data[sample]) {
      const row = tbody.insertRow(-1);
      row.insertCell(0).innerHTML = sample;
      row.insertCell(1).innerHTML = hashing_barcode;
      row.insertCell(2).innerHTML = Intl.NumberFormat("en-US").format(data[sample][hashing_barcode].n_correct + data[sample][hashing_barcode].n_corrected + data[sample][hashing_barcode].n_correct_upstream);
      row.insertCell(3).innerHTML = Intl.NumberFormat("en-US").format(data[sample][hashing_barcode].n_correct);
      row.insertCell(4).innerHTML = Intl.NumberFormat("en-US").format(data[sample][hashing_barcode].n_corrected);
      row.insertCell(5).innerHTML = Intl.NumberFormat("en-US").format(data[sample][hashing_barcode].n_correct_upstream);
    }
  }
}

document.addEventListener("DOMContentLoaded", function () {
  const binLabels = data.hashing_summary_bin_labels || [];

  generate_hashing_summary_table("sample-hashing-summary-table", data.hashing_summary || []);
  generate_hashing_summary_table("sample-hashing-summary-filt-table", data.hashing_summary_filt || []);
  generate_hashing_bins_heatmap(data.hashing_summary_bins || [], binLabels);
  generate_hashing_table(data.hashing);

  $("#sample-hashing-summary-table").tablesorter({
    sortList: [[1, 0]],
  });
  $("#sample-hashing-summary-filt-table").tablesorter({
    sortList: [[1, 0]],
  });

  $("#sample-hashing-table").tablesorter({
    sortList: [[0, 0]],
  });

  // Plotly can be initialized while the Hashing tab is hidden; force a resize when shown.
  const hashingTabToggle = document.querySelector('a[href="#tabs-hashing"]');
  if (hashingTabToggle && typeof Plotly !== "undefined") {
    hashingTabToggle.addEventListener("shown.bs.tab", function () {
      const heatmap = document.getElementById("chart-hashing-bins-heatmap");
      if (heatmap) {
        Plotly.Plots.resize(heatmap);
      }
    });
  }

});

//--------------------------------------------
// Table - Benchmarks summary (like Hashing table).
//--------------------------------------------

function generate_benchmarks_table(benchmarks) {
  const table = document.getElementById("benchmarks-table");
  if (!table) return;

  // Give it the same behavior hooks as hashing
  table.className = "table card-table table-vcenter text-nowrap datatable tablesorter";

  // Base skeleton (same style approach as hashing)
  table.innerHTML = `
      <thead>
        <tr>
          <th>Job</th>
          <th>Time<br><sub>(seconds)</sub></th>
          <th>Time<br><sub>(h:m:s)</sub></th>
          <th>Max RAM<br><sub>(gb)</sub></th>
          <th>IO Read<br><sub>(gb)</sub></th>
          <th>IO Write<br><sub>(gb)</sub></th>
          <th>Mean Load<br><sub>(cpus)</sub></th>
        </tr>
      </thead>
      <tbody class="table-tbody"></tbody>
    `;

  const tbody = table.querySelector("tbody");

  // No data case
  if (!Array.isArray(benchmarks) || benchmarks.length === 0) {
    const row = tbody.insertRow(-1);
    const cell = row.insertCell(0);
    cell.colSpan = 7;
    cell.style.textAlign = "center";
    cell.innerHTML = "<b>No benchmark data available</b>";
    return;
  }

  const fmt_number = (v, digits = 1,  divide_by = 1) => {
    if (v == null) return "N/A";
    if (typeof v === "number") {
      if (Number.isNaN(v)) return "N/A";
      return (v/divide_by).toFixed(digits);
    }
    return String(v);
  };

  // Fill rows
  for (const b of benchmarks) {
    const row = tbody.insertRow(-1);

    // Job name
    row.insertCell(-1).textContent = b.job ?? "N/A";

    // Time (seconds)
    row.insertCell(-1).textContent = fmt_number(b.s, 0);

    // Time (h:m:s)
    row.insertCell(-1).textContent = b["h:m:s"] ?? "N/A";

    // MB values that should be converted to GB for display
    const fields = ["max_rss","io_in","io_out"];
    for (const key of fields) {
      row.insertCell(-1).textContent = fmt_number(b[key], 2, 1024);
    }

    // Load (convert percentage to cpus)
    row.insertCell(-1).textContent = fmt_number(b.mean_load, 1, 100);
  }

}

document.addEventListener("DOMContentLoaded", function () {
  generate_benchmarks_table(data.benchmarks);

  $("#benchmarks-table").tablesorter({
    sortList: [[1, 1]],  // default: seconds descending
  });
});



//--------------------------------------------
// Chart - Sankey diagram of barcodes.
//--------------------------------------------

function build_uncorrectables_sankey_links(uncorrectables) {
  if(!Array.isArray(uncorrectables)) return [];
  // keys and user-visible labels for each entry in the "(True, False, ...)" tuple.
  // idx tracks the position in the original 4-element boolean tuple.
  const dims = [
    { key: "p5", label: "p5", idx: 0 },
    { key: "p7", label: "p7", idx: 1 },
    { key: "ligation", label: "ligation", idx: 2 },
    { key: "rt", label: "RT", idx: 3 },
  ].filter(d => !(d.key === "p5" && data.ignore_p5) && !(d.key === "p7" && data.ignore_p7));

  // Helper to parse tuple strings like "(True, False, True, True)" into [true, false, true, true]
  const parseTuple = (tupleStr) =>
    String(tupleStr)
      .replace(/[()]/g, "")
      .split(",")
      .map((s) => s.trim() === "True");

  const links = [];
  for (const item of uncorrectables) {
    const value = Number(item?.value ?? 0);
    const bools = parseTuple(item?.source);
    if (bools.length !== 4) continue;

    // Unique group ID for this source
    const groupId = String(item.source).replace(/\s+/g, "");

    let from = "Total Bad reads";
    for (let i = 0; i < dims.length; i++) {
      const to = `${bools[dims[i].idx] ? "Good" : "Bad"} ${dims[i].label}`;
      links.push({ from, to, value, id: `${groupId}-${i}` });
      from = to;
    }
    links.push({ from, to: " ", value, id: `${groupId}-end` });
  }

  return links;
}
const data_sankey_barcodes = build_uncorrectables_sankey_links(data.uncorrectables_sankey);

document.addEventListener("DOMContentLoaded", function () {
  document.getElementById("chart-sankey").innerHTML = ''
  am5.ready(function () {
    // Initialize root element.
    var root = am5.Root.new("chart-sankey");
    root.setThemes([am5themes_Animated.new(root)]);

    // Create series.
    var series = root.container.children.push(
      am5flow.Sankey.new(root, {
        sourceIdField: "from",
        targetIdField: "to",
        valueField: "value",
        paddingRight: 75,
        nodePadding: 20,
        idField: "id",
      })
    );

    // Add gradient to links.
    series.links.template.setAll({
      fillStyle: "gradient",
      stroke: am5.color(0x000000),
    });

    // Theme of nodes.
    series.nodes.rectangles.template.setAll({
      fillOpacity: 0.5,
      stroke: am5.color(0x000000),
      strokeWidth: 1,
      cornerRadiusTL: 4,
      cornerRadiusTR: 4,
      cornerRadiusBL: 4,
      cornerRadiusBR: 4,
    });

    // Add the values to the nodes.
    series.nodes.labels.template.setAll({
      text: "[bold]{name}[/]\n({sumOutgoing})",
    });

    // Highlight all links with the same id beginning.
    series.links.template.events.on("pointerover", function (event) {
      var dataItem = event.target.dataItem;
      var parentId = dataItem.get("id").split("-")[0];

      am5.array.each(series.dataItems, function (dataItem) {
        if (dataItem.get("id").startsWith(parentId)) {
          dataItem.get("link").hover();
        }
      });
    });

    series.links.template.events.on("pointerout", function (event) {
      am5.array.each(series.dataItems, function (dataItem) {
        dataItem.get("link").unhover();
      });
    });

    // Set data.
    series.data.setAll(data_sankey_barcodes);

    // Animation.
    series.appear(1000, 100);
  });
});
