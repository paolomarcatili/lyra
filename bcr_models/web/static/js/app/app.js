angular.module("IgSf", ['ngSanitize'])


.config(function(){


})

.controller("MainCtrl", function($scope, $http) {

    $scope.tcrForm = {};
    $scope.showForm = true;
    $scope.showOutput = false;
    $scope.outputActive = "templates"

    $scope.projectTitle = "project name";

    $scope.alphaChain = "";
    $scope.betaChain = "";

    //Heavy/Light
    //$scope.alphaChain = "EVQLVESGPGLVQPGKSLRLSCVASGFTFSGYGMHWVRQAPGKGLEWIALIIYDESNKYYADSVKGRFTISRDNSKNTLYLQMSSLRAEDTAVFYCAKVKFYDPTAPNDYWGQGTLVTVSS";
    //$scope.betaChain = "QSVLTQPPSASGTPGQRISISCSGTSSNVENNYVYWYQHLPGTAPKLLIYRNDHRSSGIPDRFSASKSGTSASLAISGLRPEDEGDYYCAAWDDSRGGPDWVFGGGTKLTVLAQP";

    $scope.PDB = ""

    var defaultInputStatus = 'Input alpha and beta chain protein sequences..';
    $scope.inputStatus = defaultInputStatus

    $scope.submitAuto = function(item, event) {
        if ($scope.alphaChain && $scope.betaChain) {
            var dataObject = {
                alphaChain: $scope.alphaChain,
                betaChain: $scope.betaChain
            }

            var responsePromise = $http.post("/model", dataObject, {});
            responsePromise.success(function(dataFromServer, status, headers, config) {
                $scope.JSON = dataFromServer;
                $scope.showForm = false;
                $scope.showOutput = true;
            });
            responsePromise.error(function(data, status, headers, config) {
                alert("Submitting form failed!");
            });
        }
     }

    $scope.download = function() {
        saveTextAs($scope.JSON.pdb, "structure.pdb");
     }

    $scope.submitManual = function(item, event) {
        alert('Not implemented..');
    }

    $scope.loadExample = function() {
        $scope.alphaChain = "DSVTQMQGQVTLSENDFLFINCTYSTTGYPTLFWYVQYSGEGPQLLLQVTTANNKGSSRGFEATYDKGTTSFHLQKTSVQEIDSAVYYCAANSGTYQRFGTGTKLQVVP";
        $scope.betaChain = "AVTQSPRNKVAVTGGKVTLSCNQTNNHNNMYWYRQDTGHGLRLIHYSYGAGSTEKGDIPDGYKASRPSQENFSLILELATPSQTSVYFCASGDFWGDTLYFGAGTRLSVL";
    }


    $scope.invalidChains = function() {
        if ($scope.alphaChain && $scope.betaChain){
            $scope.inputStatus = 'Ready for modelling..';
            return false;
        } else {
            $scope.inputStatus = defaultInputStatus;
            return true;
        }
     }

});