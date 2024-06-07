function showTooltip(evt, text) {
  let tooltip = document.getElementById("tooltip");
  tooltip.innerHTML = text;
  tooltip.style.display = "block";
  tooltip.style.left = evt.pageX + 10 + 'px';
  tooltip.style.top = evt.pageY + 10 + 'px';
}

function hideTooltip() {
  var tooltip = document.getElementById("tooltip");
  tooltip.style.display = "none";
}

// var myicon = document.getElementById("myicon");
// var mypopup = document.getElementById("mypopup");

// myicon.addEventListener("mouseover", showPopup);
// myicon.addEventListener("mouseout", hidePopup);

// function showPopup(evt) {
//   var iconPos = myicon.getBoundingClientRect();
//   mypopup.style.left = (iconPos.right + 20) + "px";
//   mypopup.style.top = (window.scrollY + iconPos.top - 60) + "px";
//   mypopup.style.display = "block";
// }

// function hidePopup(evt) {
//   mypopup.style.display = "none";
// }